"""
CLI entry point for PubMed cost-ratio literature search.

Usage:
    python -m scripts.pubmed_search.main                            # all stages (no fulltext)
    python -m scripts.pubmed_search.main --stage search             # search PubMed only
    python -m scripts.pubmed_search.main --stage fulltext           # download PDFs via Sci-Hub
    python -m scripts.pubmed_search.main --stage fulltext_extract   # extract text from PDFs
    python -m scripts.pubmed_search.main --stage extract            # extract costs (uses fulltext if available)
    python -m scripts.pubmed_search.main --stage output             # build QMD
    python -m scripts.pubmed_search.main --stage all_with_fulltext  # full pipeline including Sci-Hub
    python -m scripts.pubmed_search.main --diseases cholera dengue  # subset of diseases
    python -m scripts.pubmed_search.main --max-pdfs 10              # limit PDFs per disease
"""

import argparse
import sys

from .config import DISEASES, FULLTEXT_MAP_PATH
from .search_runner import search_all_diseases
from .extract_costs import extract_all
from .build_output import build_qmd


def main():
    parser = argparse.ArgumentParser(
        description="PubMed literature search for vaccine cost ratios (R = C_I/C_V)"
    )
    parser.add_argument(
        "--stage",
        choices=[
            "search", "fulltext", "fulltext_extract",
            "extract", "output", "all", "all_with_fulltext",
        ],
        default="all",
        help="Pipeline stage to run (default: all = search+extract+output)",
    )
    parser.add_argument(
        "--diseases",
        nargs="+",
        default=None,
        help=f"Disease keys to process. Available: {', '.join(DISEASES.keys())}",
    )
    parser.add_argument(
        "--max-pdfs",
        type=int,
        default=30,
        help="Max PDFs to download per disease (default: 30)",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=5.0,
        help="Seconds between Sci-Hub downloads (default: 5)",
    )
    args = parser.parse_args()

    disease_keys = args.diseases
    if disease_keys:
        unknown = [k for k in disease_keys if k not in DISEASES]
        if unknown:
            print(f"ERROR: Unknown disease keys: {unknown}")
            print(f"Available: {', '.join(DISEASES.keys())}")
            sys.exit(1)

    if args.stage == "all":
        stages = ["search", "extract", "output"]
    elif args.stage == "all_with_fulltext":
        stages = ["search", "fulltext", "fulltext_extract", "extract", "output"]
    else:
        stages = [args.stage]

    if "search" in stages:
        print("=" * 60)
        print("STAGE 1: PubMed Search")
        print("=" * 60)
        search_all_diseases(disease_keys=disease_keys)
        print()

    if "fulltext" in stages:
        print("=" * 60)
        print("STAGE 1.5: Download Full-Text PDFs (Sci-Hub)")
        print("=" * 60)
        from .fulltext_fetcher import download_fulltext

        download_fulltext(
            disease_keys=disease_keys,
            max_per_disease=args.max_pdfs,
            delay=args.delay,
        )
        print()

    if "fulltext_extract" in stages:
        print("=" * 60)
        print("STAGE 1.6: Extract Text from PDFs")
        print("=" * 60)
        from .fulltext_extractor import extract_fulltext_for_cache
        from .config import RAW_CACHE_PATH

        extract_fulltext_for_cache(
            cache_path=RAW_CACHE_PATH,
            output_path=FULLTEXT_MAP_PATH,
        )
        print()

    if "extract" in stages:
        print("=" * 60)
        print("STAGE 2: Cost Data Extraction")
        print("=" * 60)
        # Use fulltext map if it exists
        import os

        ft_path = FULLTEXT_MAP_PATH if os.path.exists(FULLTEXT_MAP_PATH) else None
        extract_all(disease_keys=disease_keys, fulltext_path=ft_path)
        print()

    if "output" in stages:
        print("=" * 60)
        print("STAGE 3: Build Quarto Document")
        print("=" * 60)
        build_qmd()
        print()

    print("Done.")


if __name__ == "__main__":
    main()
