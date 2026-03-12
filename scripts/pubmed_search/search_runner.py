"""
Stage 1: Search PubMed for cost/economic studies by disease and cache results.
"""

import json
import os
from pathlib import Path

from .config import DISEASES, COST_TERMS, MAX_RESULTS_PER_DISEASE, RAW_CACHE_PATH
from .pubmed_api import esearch, efetch_abstracts


def build_query(disease_key: str) -> str:
    """
    Build a PubMed query combining disease synonyms with economic cost terms.

    The query structure is:
        (synonym1 OR synonym2 OR ...) AND (cost_term1 OR cost_term2 OR ...)
    """
    disease = DISEASES[disease_key]
    disease_part = " OR ".join(f'"{s}"' for s in disease["synonyms"])
    cost_part = " OR ".join(f'"{t}"' for t in COST_TERMS)
    return f"({disease_part}) AND ({cost_part})"


def search_disease(disease_key: str) -> dict:
    """
    Search PubMed for one disease and return structured results.

    Returns
    -------
    dict
        Keys: disease_key, disease_name, query, n_results, articles.
    """
    query = build_query(disease_key)
    print(f"  Searching: {disease_key} ...", end=" ", flush=True)

    pmids = esearch(query, retmax=MAX_RESULTS_PER_DISEASE)
    print(f"{len(pmids)} hits.", end=" ", flush=True)

    if not pmids:
        print("No abstracts to fetch.")
        return {
            "disease_key": disease_key,
            "disease_name": DISEASES[disease_key]["name"],
            "query": query,
            "n_results": 0,
            "articles": [],
        }

    articles = efetch_abstracts(pmids)
    # Keep only articles that actually have an abstract
    articles = [a for a in articles if a["abstract"].strip()]
    print(f"{len(articles)} with abstracts.")

    return {
        "disease_key": disease_key,
        "disease_name": DISEASES[disease_key]["name"],
        "query": query,
        "n_results": len(articles),
        "articles": articles,
    }


def search_all_diseases(
    cache_path: str = RAW_CACHE_PATH,
    disease_keys: list[str] | None = None,
) -> dict:
    """
    Search PubMed for all (or selected) diseases and save JSON cache.

    Parameters
    ----------
    cache_path : str
        Where to write the JSON cache (relative to project root).
    disease_keys : list[str] or None
        Subset of disease keys to search. None = all 14.

    Returns
    -------
    dict
        Full results dict keyed by disease_key.
    """
    keys = disease_keys or list(DISEASES.keys())

    # Load existing cache if present (to allow incremental runs)
    results = {}
    if os.path.exists(cache_path):
        with open(cache_path, "r", encoding="utf-8") as f:
            results = json.load(f)
        print(f"Loaded existing cache with {len(results)} diseases.")

    for key in keys:
        if key not in DISEASES:
            print(f"  WARNING: unknown disease key '{key}', skipping.")
            continue
        result = search_disease(key)
        results[key] = result

    # Save
    Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"Cache saved to {cache_path}")

    return results
