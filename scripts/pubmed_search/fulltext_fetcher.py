"""
Stage 1.5: Download full-text PDFs via Sci-Hub for cached PubMed articles.

Uses DOI (preferred) or PMID as identifier. Adds a configurable delay
between downloads to avoid aggressive rate limiting. Stores PDFs in
data/fulltext_pdfs/{disease_key}/ and updates the JSON cache with
download status.
"""

import json
import os
import time
from pathlib import Path

from .config import RAW_CACHE_PATH
from .scihub_client import SciHubClient

PDF_DIR = "data/fulltext_pdfs"
DOWNLOAD_DELAY = 5  # seconds between downloads (be polite)
MAX_PER_DISEASE = 30  # limit downloads per disease to avoid blocks


def _get_doi_from_article(article: dict) -> str | None:
    """Try to extract DOI from article metadata or abstract text."""
    # Check if stored in article dict (we add it during efetch parsing)
    doi = article.get("doi")
    if doi:
        return doi
    # Some abstracts mention DOI
    import re
    match = re.search(r"10\.\d{4,9}/[^\s]+", article.get("abstract", ""))
    if match:
        return match.group(0).rstrip(".,;)")
    return None


def download_fulltext(
    cache_path: str = RAW_CACHE_PATH,
    pdf_dir: str = PDF_DIR,
    disease_keys: list[str] | None = None,
    max_per_disease: int = MAX_PER_DISEASE,
    delay: float = DOWNLOAD_DELAY,
    priority: str = "cost_relevant",
) -> dict:
    """
    Download full-text PDFs for cached articles.

    Parameters
    ----------
    cache_path : str
        Path to PubMed raw results JSON cache.
    pdf_dir : str
        Directory to store downloaded PDFs.
    disease_keys : list[str] or None
        Subset of diseases. None = all in cache.
    max_per_disease : int
        Maximum PDFs to download per disease.
    delay : float
        Seconds between download attempts.
    priority : str
        "cost_relevant" — prioritise articles mentioning cost/economic terms.
        "all" — download in order.

    Returns
    -------
    dict
        Summary: {disease_key: {attempted, success, failed, skipped}}.
    """
    with open(cache_path, "r", encoding="utf-8") as f:
        cache = json.load(f)

    keys = disease_keys or list(cache.keys())
    client = SciHubClient(timeout=30)
    summary = {}

    for disease_key in keys:
        if disease_key not in cache:
            print(f"  WARNING: '{disease_key}' not in cache, skipping.")
            continue

        disease_data = cache[disease_key]
        articles = disease_data.get("articles", [])
        disease_dir = os.path.join(pdf_dir, disease_key)
        os.makedirs(disease_dir, exist_ok=True)

        # Prioritise articles that mention cost terms
        if priority == "cost_relevant":
            import re
            cost_re = re.compile(
                r"cost|economic|expenditure|burden|DALY|QALY|"
                r"cost-effective|budget|financial",
                re.IGNORECASE,
            )
            articles = sorted(
                articles,
                key=lambda a: bool(cost_re.search(a.get("abstract", ""))),
                reverse=True,
            )

        stats = {"attempted": 0, "success": 0, "failed": 0, "skipped": 0}

        for article in articles:
            if stats["success"] >= max_per_disease:
                break

            pmid = article.get("pmid", "")
            pdf_path = os.path.join(disease_dir, f"{pmid}.pdf")

            # Skip if already downloaded
            if os.path.exists(pdf_path) and os.path.getsize(pdf_path) > 1000:
                stats["skipped"] += 1
                article["fulltext_pdf"] = pdf_path
                continue

            # DOI is required — Sci-Hub mirrors (especially .kr) need DOI
            doi = _get_doi_from_article(article)
            if not doi:
                stats["failed"] += 1
                continue

            stats["attempted"] += 1
            print(
                f"  [{disease_key}] PMID {pmid} (DOI: {doi}) ... ",
                end="",
                flush=True,
            )

            result = client.fetch(doi)

            if "err" not in result:
                with open(pdf_path, "wb") as f:
                    f.write(result["pdf"])
                article["fulltext_pdf"] = pdf_path
                stats["success"] += 1
                print("OK")
            else:
                stats["failed"] += 1
                print(f"FAILED: {result['err'][:80]}")

            time.sleep(delay)

        summary[disease_key] = stats
        print(
            f"  {disease_key}: {stats['success']} downloaded, "
            f"{stats['failed']} failed, {stats['skipped']} skipped"
        )

    # Update cache with fulltext_pdf paths
    with open(cache_path, "w", encoding="utf-8") as f:
        json.dump(cache, f, indent=2, ensure_ascii=False)
    print(f"Cache updated with PDF paths at {cache_path}")

    return summary
