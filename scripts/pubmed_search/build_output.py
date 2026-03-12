"""
Stage 3: Build Quarto document from extracted cost data.
"""

import json
from pathlib import Path
from collections import defaultdict

from .config import DISEASES, EXTRACTION_PATH, QMD_OUTPUT_PATH


def _summarise_disease(extractions: list[dict]) -> dict:
    """Summarise extractions for a single disease."""
    n_articles = len(extractions)
    years = [e["year"] for e in extractions if e["year"]]
    journals = [e["journal"] for e in extractions if e["journal"]]
    countries = []
    for e in extractions:
        countries.extend(e.get("countries", []))

    study_types = defaultdict(int)
    for e in extractions:
        for st in e.get("study_types", []):
            study_types[st] += 1

    # Collect cost-per-case mentions
    cost_per_case_mentions = []
    vaccine_cost_mentions = []
    for e in extractions:
        for m in e.get("money_mentions", []):
            if "cost_per_case" in m.get("components", []):
                cost_per_case_mentions.append(
                    {
                        "amount": m["amount_usd"],
                        "raw": m["amount_raw"],
                        "context": m["context"],
                        "pmid": e["pmid"],
                        "year": e["year"],
                    }
                )
            if "vaccine_cost" in m.get("components", []) or "delivery_cost" in m.get(
                "components", []
            ):
                vaccine_cost_mentions.append(
                    {
                        "amount": m["amount_usd"],
                        "raw": m["amount_raw"],
                        "context": m["context"],
                        "pmid": e["pmid"],
                        "year": e["year"],
                    }
                )

    # Unique countries, deduped
    unique_countries = list(dict.fromkeys(countries))

    return {
        "n_articles": n_articles,
        "year_range": f"{min(years, default='?')}–{max(years, default='?')}",
        "study_types": dict(study_types),
        "countries": unique_countries[:10],
        "cost_per_case": cost_per_case_mentions,
        "vaccine_cost": vaccine_cost_mentions,
    }


def build_qmd(
    extraction_path: str = EXTRACTION_PATH,
    output_path: str = QMD_OUTPUT_PATH,
) -> str:
    """
    Generate a Quarto document summarising PubMed cost-ratio literature.

    Returns the path to the generated file.
    """
    with open(extraction_path, "r", encoding="utf-8") as f:
        extractions = json.load(f)

    # Group by disease
    by_disease = defaultdict(list)
    for e in extractions:
        by_disease[e["disease_key"]].append(e)

    lines = []
    lines.append("---")
    lines.append('title: "Cost Ratio Literature Review: R = C_I / C_V"')
    lines.append('subtitle: "PubMed-based survey across VIMC diseases"')
    lines.append("format:")
    lines.append("  html:")
    lines.append("    toc: true")
    lines.append("    toc-depth: 3")
    lines.append("    code-fold: true")
    lines.append("    embed-resources: true")
    lines.append("---")
    lines.append("")
    lines.append("## Overview")
    lines.append("")
    lines.append(
        "This document summarises cost data extracted from PubMed abstracts "
        "for 14 vaccine-preventable diseases listed by the Vaccine Impact "
        "Modelling Consortium (VIMC). The goal is to estimate plausible ranges "
        "of the cost ratio $R = C_I / C_V$ (infection cost per case divided by "
        "vaccination cost per person) used in the reactive vs. pre-emptive "
        "vaccination framework."
    )
    lines.append("")
    lines.append(
        "::: {.callout-note}\n"
        "Cost data were extracted from **abstracts only** using regex heuristics. "
        "Many cost figures appear only in full-text articles. This survey is a "
        "starting point for manual curation, not a definitive extraction.\n"
        ":::"
    )
    lines.append("")

    # Summary table
    lines.append("## Summary by Disease")
    lines.append("")
    lines.append(
        "| Disease | Articles | Year range | Study types | "
        "Cost-per-case mentions | Vaccine cost mentions |"
    )
    lines.append("|---------|----------|------------|-------------|"
                 "----------------------|----------------------|")

    for key in DISEASES:
        if key not in by_disease:
            lines.append(
                f"| {DISEASES[key]['name']} | 0 | — | — | 0 | 0 |"
            )
            continue
        summary = _summarise_disease(by_disease[key])
        st_str = ", ".join(
            f"{k}({v})" for k, v in summary["study_types"].items()
        )
        lines.append(
            f"| {DISEASES[key]['name']} | {summary['n_articles']} | "
            f"{summary['year_range']} | {st_str} | "
            f"{len(summary['cost_per_case'])} | "
            f"{len(summary['vaccine_cost'])} |"
        )

    lines.append("")

    # Per-disease sections
    for key in DISEASES:
        if key not in by_disease:
            continue
        disease_name = DISEASES[key]["name"]
        exts = by_disease[key]
        summary = _summarise_disease(exts)

        lines.append(f"## {disease_name}")
        lines.append("")
        lines.append(
            f"**{summary['n_articles']} articles** with cost data "
            f"({summary['year_range']}). "
            f"Countries: {', '.join(summary['countries'][:5]) or 'not specified'}."
        )
        lines.append("")

        # Cost per case table
        if summary["cost_per_case"]:
            lines.append("### Cost per case mentions")
            lines.append("")
            lines.append("| Amount | PMID | Year | Context (excerpt) |")
            lines.append("|--------|------|------|-------------------|")
            for m in summary["cost_per_case"][:15]:
                ctx = m["context"][:120].replace("|", "\\|").replace("\n", " ")
                lines.append(
                    f"| {m['raw']} | [{m['pmid']}]"
                    f"(https://pubmed.ncbi.nlm.nih.gov/{m['pmid']}/) | "
                    f"{m['year']} | {ctx} |"
                )
            lines.append("")

        # Vaccine cost table
        if summary["vaccine_cost"]:
            lines.append("### Vaccine cost mentions")
            lines.append("")
            lines.append("| Amount | PMID | Year | Context (excerpt) |")
            lines.append("|--------|------|------|-------------------|")
            for m in summary["vaccine_cost"][:10]:
                ctx = m["context"][:120].replace("|", "\\|").replace("\n", " ")
                lines.append(
                    f"| {m['raw']} | [{m['pmid']}]"
                    f"(https://pubmed.ncbi.nlm.nih.gov/{m['pmid']}/) | "
                    f"{m['year']} | {ctx} |"
                )
            lines.append("")

        # Top articles for manual review
        lines.append("### Key articles for review")
        lines.append("")
        # Prioritise articles with both cost-per-case and vaccine cost mentions
        scored = []
        for e in exts:
            score = e["n_money_mentions"] + (3 if e["has_cost_per_case"] else 0) + (
                2 if e["has_vaccine_cost"] else 0
            )
            scored.append((score, e))
        scored.sort(key=lambda x: -x[0])

        for _, e in scored[:5]:
            pmid = e["pmid"]
            lines.append(
                f"- [{e['title']}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/) "
                f"({e['year']}, {e['journal']}). "
                f"Types: {', '.join(e['study_types'])}. "
                f"Money mentions: {e['n_money_mentions']}."
            )
        lines.append("")

    # Methodology note
    lines.append("## Methodology")
    lines.append("")
    lines.append(
        "Articles were retrieved via NCBI E-utilities (PubMed) using disease-specific "
        "search terms combined with economic/cost keywords. Cost data were extracted "
        "from abstracts using regex pattern matching for USD amounts. Study types were "
        "classified by keyword matching (CEA, COI, CBA, BIA). Countries/regions were "
        "identified by name matching."
    )
    lines.append("")
    lines.append(
        "**Limitations**: (1) Only abstracts searched — many cost figures are in "
        "full text only. (2) Currency conversion and inflation adjustment not applied. "
        "(3) Regex extraction is approximate — manual verification needed. "
        "(4) Non-English articles excluded."
    )

    content = "\n".join(lines)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)

    print(f"Quarto document saved to {output_path}")
    return output_path
