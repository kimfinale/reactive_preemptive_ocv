"""
Vendored and patched Sci-Hub client.

Based on https://github.com/zaytoun/scihub.py (MIT license, @zaytoun).
Patches:
  - Hardcoded fallback Sci-Hub URLs (the upstream URL discovery is broken).
  - Reduced retry noise, added timeout controls.
  - Simplified for our use case (fetch by DOI/PMID only, no Google Scholar).
"""

import hashlib
import logging
import re

import requests
import urllib3
from bs4 import BeautifulSoup

urllib3.disable_warnings()
logger = logging.getLogger("scihub_client")

HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/120.0.0.0 Safari/537.36"
    )
}

# Fallback Sci-Hub mirrors (update as needed — checked March 2026)
SCIHUB_URLS = [
    "https://sci-hub.kr",
    "https://sci-hub.st",
    "https://sci-hub.se",
    "https://sci-hub.ru",
    "https://sci-hub.red",
    "https://sci-hub.box",
    "https://sci-hub.ee",
]


class CaptchaNeedException(Exception):
    pass


class SciHubClient:
    """
    Minimal Sci-Hub client for fetching papers by DOI or PMID.
    """

    def __init__(self, base_urls: list[str] | None = None, timeout: int = 30):
        self.sess = requests.Session()
        self.sess.headers.update(HEADERS)
        self._all_urls = list(base_urls or SCIHUB_URLS)
        self._mirror_idx = 0
        self.base_url = self._all_urls[0]
        self.timeout = timeout

    def _change_base_url(self):
        self._mirror_idx += 1
        if self._mirror_idx >= len(self._all_urls):
            raise RuntimeError("Ran out of Sci-Hub mirrors to try.")
        self.base_url = self._all_urls[self._mirror_idx]
        logger.info("Switching to %s", self.base_url)

    def _reset_mirrors(self):
        """Reset mirror index to start from the first mirror."""
        self._mirror_idx = 0
        self.base_url = self._all_urls[0]

    def fetch(self, identifier: str) -> dict:
        """
        Fetch a paper by DOI, PMID, or direct URL.

        Returns dict with keys:
          - 'pdf': bytes (raw PDF content)
          - 'url': str (resolved URL)
          - 'name': str (generated filename)
        Or on failure:
          - 'err': str (error message)
        """
        self._reset_mirrors()
        last_err = ""
        n_mirrors = len(self._all_urls)
        for attempt in range(n_mirrors):
            try:
                return self._try_fetch(identifier)
            except CaptchaNeedException as e:
                last_err = str(e)
                logger.warning("CAPTCHA on %s, trying next mirror", self.base_url)
            except requests.exceptions.HTTPError as e:
                last_err = str(e)
                logger.warning("HTTP %s on %s", e.response.status_code if e.response else "?", self.base_url)
            except requests.exceptions.ConnectionError as e:
                last_err = str(e)
                logger.warning("Connection failed on %s", self.base_url)
            except requests.exceptions.RequestException as e:
                last_err = str(e)
                logger.warning("Request error on %s: %s", self.base_url, e)

            # Try next mirror if available
            if attempt < n_mirrors - 1:
                try:
                    self._change_base_url()
                except RuntimeError:
                    break

        return {"err": f"All mirrors failed for {identifier}: {last_err}"}

    def _try_fetch(self, identifier: str) -> dict:
        """Single attempt to fetch from current base_url."""
        # If it's already a direct PDF URL, download directly
        if identifier.startswith("http") and identifier.endswith(".pdf"):
            res = self.sess.get(identifier, verify=False, timeout=self.timeout)
            res.raise_for_status()
            return {
                "pdf": res.content,
                "url": identifier,
                "name": self._generate_name(res),
            }

        # Otherwise, go through Sci-Hub
        url = f"{self.base_url}/{identifier}"
        res = self.sess.get(url, verify=False, timeout=self.timeout)
        res.raise_for_status()  # raises HTTPError → caught by fetch() to try next mirror

        # Parse the page for the PDF link
        soup = BeautifulSoup(res.content, "html.parser")

        pdf_url = None

        # Method 1: citation_pdf_url meta tag (used by sci-hub.kr)
        meta = soup.find("meta", attrs={"name": "citation_pdf_url"})
        if meta and meta.get("content"):
            pdf_url = meta["content"]

        # Method 2: iframe or embed tag
        if not pdf_url:
            for tag_name in ("iframe", "embed"):
                tag = soup.find(tag_name)
                if tag and tag.get("src"):
                    pdf_url = tag["src"]
                    break

        # Method 3: direct PDF link
        if not pdf_url:
            for a in soup.find_all("a", href=True):
                if a["href"].endswith(".pdf"):
                    pdf_url = a["href"]
                    break
            # Check onclick with location.href
            btn = soup.find("button", onclick=True)
            if btn and ".pdf" in btn.get("onclick", ""):
                match = re.search(r"location\.href\s*=\s*'([^']+)'", btn["onclick"])
                if match:
                    pdf_url = match.group(1)

        # If no PDF link found, check if it's a CAPTCHA or just not found
        if not pdf_url:
            if "captcha" in res.text.lower() and "notfound" not in res.text.lower():
                raise CaptchaNeedException(f"CAPTCHA for {identifier}")
            return {"err": f"Could not find PDF link for {identifier} on {self.base_url}"}

        # Normalize URL
        if pdf_url.startswith("//"):
            pdf_url = "https:" + pdf_url
        elif not pdf_url.startswith("http"):
            pdf_url = self.base_url + pdf_url

        # Download the actual PDF
        pdf_res = self.sess.get(pdf_url, verify=False, timeout=self.timeout)
        pdf_res.raise_for_status()

        if "application/pdf" not in pdf_res.headers.get("Content-Type", ""):
            raise CaptchaNeedException(
                f"Got non-PDF content type for {identifier}: "
                f"{pdf_res.headers.get('Content-Type')}"
            )

        return {
            "pdf": pdf_res.content,
            "url": pdf_url,
            "name": self._generate_name(pdf_res),
        }

    @staticmethod
    def _generate_name(res: requests.Response) -> str:
        name = res.url.split("/")[-1]
        name = re.sub(r"#view=(.+)", "", name)
        pdf_hash = hashlib.md5(res.content).hexdigest()[:8]
        return f"{pdf_hash}-{name[-30:]}"

    def download(
        self,
        identifier: str,
        destination: str = "",
        filename: str | None = None,
    ) -> dict:
        """Fetch and save PDF to disk."""
        import os

        result = self.fetch(identifier)
        if "err" not in result:
            path = os.path.join(destination, filename or result["name"])
            os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
            with open(path, "wb") as f:
                f.write(result["pdf"])
            result["path"] = path
        return result
