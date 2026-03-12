"""
Extract text from downloaded full-text PDFs using PyMuPDF (fitz).
"""

import os

try:
    import fitz  # PyMuPDF
except ImportError:
    fitz = None


def extract_text_from_pdf(pdf_path: str, max_pages: int = 50) -> str:
    """
    Extract text from a PDF file.

    Parameters
    ----------
    pdf_path : str
        Path to PDF file.
    max_pages : int
        Maximum number of pages to extract (to cap very long documents).

    Returns
    -------
    str
        Extracted text, or empty string on failure.
    """
    if fitz is None:
        raise ImportError(
            "PyMuPDF is required for PDF text extraction. "
            "Install with: pip install pymupdf"
        )

    if not os.path.exists(pdf_path) or os.path.getsize(pdf_path) < 1000:
        return ""

    try:
        doc = fitz.open(pdf_path)
        pages = []
        for i, page in enumerate(doc):
            if i >= max_pages:
                break
            pages.append(page.get_text())
        doc.close()
        return "\n".join(pages)
    except Exception as e:
        print(f"  WARNING: Failed to extract text from {pdf_path}: {e}")
        return ""


def extract_fulltext_for_cache(
    cache_path: str,
    output_path: str | None = None,
) -> dict:
    """
    Extract text from all PDFs referenced in the cache and store
    the text alongside each article entry.

    Parameters
    ----------
    cache_path : str
        Path to PubMed raw results JSON cache (with fulltext_pdf paths).
    output_path : str or None
        If provided, save a separate JSON with {pmid: fulltext} mapping.

    Returns
    -------
    dict
        {pmid: fulltext_text} for all successfully extracted articles.
    """
    import json

    with open(cache_path, "r", encoding="utf-8") as f:
        cache = json.load(f)

    fulltext_map = {}
    total = 0
    success = 0

    for disease_key, disease_data in cache.items():
        for article in disease_data.get("articles", []):
            pdf_path = article.get("fulltext_pdf")
            if not pdf_path:
                continue

            total += 1
            pmid = article.get("pmid", "")
            text = extract_text_from_pdf(pdf_path)

            if text.strip():
                # Store just the length in cache to keep it manageable
                article["fulltext_length"] = len(text)
                fulltext_map[pmid] = text
                success += 1

    print(f"Extracted text from {success}/{total} PDFs.")

    # Save fulltext map separately (can be large)
    if output_path:
        from pathlib import Path

        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(fulltext_map, f, ensure_ascii=False)
        print(f"Fulltext map saved to {output_path}")

    # Update cache
    with open(cache_path, "w", encoding="utf-8") as f:
        json.dump(cache, f, indent=2, ensure_ascii=False)

    return fulltext_map
