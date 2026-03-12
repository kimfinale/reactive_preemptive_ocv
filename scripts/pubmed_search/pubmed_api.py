"""
NCBI E-utilities wrapper for PubMed search and abstract retrieval.
Rate-limited to comply with NCBI usage policy (< 3 requests/sec).
"""

import time
import xml.etree.ElementTree as ET
from typing import Optional

import requests

from .config import (
    ESEARCH_URL,
    EFETCH_URL,
    RATE_LIMIT_DELAY,
    EFETCH_BATCH_SIZE,
)

_last_request_time = 0.0


def _rate_limit():
    """Enforce minimum delay between NCBI requests."""
    global _last_request_time
    elapsed = time.time() - _last_request_time
    if elapsed < RATE_LIMIT_DELAY:
        time.sleep(RATE_LIMIT_DELAY - elapsed)
    _last_request_time = time.time()


def esearch(query: str, retmax: int = 200) -> list[str]:
    """
    Search PubMed and return a list of PMIDs.

    Parameters
    ----------
    query : str
        PubMed search query string.
    retmax : int
        Maximum number of results to return.

    Returns
    -------
    list[str]
        List of PMID strings.
    """
    _rate_limit()
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": retmax,
        "retmode": "json",
        "sort": "relevance",
    }
    resp = requests.get(ESEARCH_URL, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    return data.get("esearchresult", {}).get("idlist", [])


def efetch_abstracts(pmids: list[str]) -> list[dict]:
    """
    Fetch article metadata and abstracts for a list of PMIDs.

    Parameters
    ----------
    pmids : list[str]
        PubMed IDs to fetch.

    Returns
    -------
    list[dict]
        Each dict has keys: pmid, title, abstract, year, journal, mesh_terms.
    """
    articles = []
    for i in range(0, len(pmids), EFETCH_BATCH_SIZE):
        batch = pmids[i : i + EFETCH_BATCH_SIZE]
        _rate_limit()
        params = {
            "db": "pubmed",
            "id": ",".join(batch),
            "rettype": "xml",
            "retmode": "xml",
        }
        resp = requests.get(EFETCH_URL, params=params, timeout=60)
        resp.raise_for_status()
        articles.extend(_parse_pubmed_xml(resp.text))
    return articles


def _parse_pubmed_xml(xml_text: str) -> list[dict]:
    """Parse PubMed XML response into structured article dicts."""
    root = ET.fromstring(xml_text)
    articles = []

    for article_elem in root.findall(".//PubmedArticle"):
        pmid_el = article_elem.find(".//PMID")
        pmid = pmid_el.text if pmid_el is not None else ""

        title_el = article_elem.find(".//ArticleTitle")
        title = _get_text(title_el)

        # Abstract may have multiple AbstractText elements (structured abstract)
        abstract_parts = []
        for ab in article_elem.findall(".//AbstractText"):
            label = ab.get("Label", "")
            text = _get_text(ab)
            if label:
                abstract_parts.append(f"{label}: {text}")
            else:
                abstract_parts.append(text)
        abstract = " ".join(abstract_parts)

        # Year
        year_el = article_elem.find(".//PubDate/Year")
        medline_year = article_elem.find(".//PubDate/MedlineDate")
        if year_el is not None:
            year = year_el.text
        elif medline_year is not None:
            year = medline_year.text[:4]
        else:
            year = ""

        # Journal
        journal_el = article_elem.find(".//Journal/Title")
        journal = journal_el.text if journal_el is not None else ""

        # DOI
        doi = ""
        for aid in article_elem.findall(".//ArticleId"):
            if aid.get("IdType") == "doi" and aid.text:
                doi = aid.text.strip()
                break
        # Fallback: check ELocationID
        if not doi:
            for eloc in article_elem.findall(".//ELocationID"):
                if eloc.get("EIdType") == "doi" and eloc.text:
                    doi = eloc.text.strip()
                    break

        # MeSH terms
        mesh_terms = []
        for mesh in article_elem.findall(".//MeshHeading/DescriptorName"):
            if mesh.text:
                mesh_terms.append(mesh.text)

        articles.append(
            {
                "pmid": pmid,
                "doi": doi,
                "title": title,
                "abstract": abstract,
                "year": year,
                "journal": journal,
                "mesh_terms": mesh_terms,
            }
        )

    return articles


def _get_text(element: Optional[ET.Element]) -> str:
    """Extract all text content from an XML element, including mixed content."""
    if element is None:
        return ""
    return "".join(element.itertext()).strip()
