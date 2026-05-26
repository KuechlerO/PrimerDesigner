"""SNP / dbSNP awareness for allelic variant primer design (Ensembl overlap API)."""

from __future__ import annotations

import logging
from typing import Any, List, Optional, Sequence, Tuple

from primer_designer_app.utils.display_utils import compute_display_bounds
from primer_designer_app.utils.ensembl_client import (
    GNOMAD_VARIANT_SET,
    EnsemblClient,
)
from primer_designer_app.utils.primer_utils import PrimerPairResult
from primer_designer_app.utils.variant_info import (
    AllelicVariantInfo,
    SequenceVariantInfo,
    VARIANT_FLANKING,
)

LOGGER = logging.getLogger(__name__)

# Common-variant threshold (gnomAD / global MAF); strictly greater than 1%.
COMMON_VARIANT_MAF_THRESHOLD = 0.01
GNOMAD_GLOBAL_POPULATIONS = ("gnomADg:ALL", "1000GENOMES:phase_3:ALL")

SNP_STATUS_NONE = "none"
SNP_STATUS_CAUTION = "caution"
SNP_STATUS_CONFLICT = "conflict"
SNP_STATUS_SKIPPED = "skipped"
SNP_STATUS_ERROR = "error"


def _intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and b_start <= a_end


def effective_maf(variation: dict) -> float | None:
    """
    Return minor allele frequency from Ensembl variation JSON (pops=1).

    Prefer the top-level MAF field; otherwise use gnomAD global or 1000 Genomes ALL.
    """
    maf = variation.get("MAF")
    if maf is not None:
        try:
            return float(maf)
        except (TypeError, ValueError):
            pass

    minor_allele = variation.get("minor_allele")
    for pop_key in GNOMAD_GLOBAL_POPULATIONS:
        allele_freqs: dict[str, float] = {}
        for entry in variation.get("populations", []):
            if entry.get("population") != pop_key:
                continue
            allele = entry.get("allele")
            freq = entry.get("frequency")
            if allele is None or freq is None:
                continue
            try:
                allele_freqs[str(allele)] = float(freq)
            except (TypeError, ValueError):
                continue
        if not allele_freqs:
            continue
        if minor_allele is not None and str(minor_allele) in allele_freqs:
            return allele_freqs[str(minor_allele)]
        return min(allele_freqs.values())

    return None


def _attach_maf_and_filter_common_variants(
    hits: List[dict],
    variation_details: dict[str, dict],
    *,
    maf_threshold: float = COMMON_VARIANT_MAF_THRESHOLD,
) -> List[dict]:
    """Keep only variants with global MAF strictly above the threshold."""
    filtered: List[dict] = []
    for hit in hits:
        var_id = hit.get("id") or ""
        details = variation_details.get(var_id)
        if not details:
            continue
        maf = effective_maf(details)
        if maf is None or maf <= maf_threshold:
            continue
        filtered.append({**hit, "maf": round(maf, 6)})
    return filtered


def get_design_region_genomic(
    var_info: AllelicVariantInfo,
    flank: int = VARIANT_FLANKING,
    primer_pairs: Optional[Sequence[PrimerPairResult]] = None,
    *,
    primer_target: Optional[Tuple[int, int]] = None,
) -> Optional[dict]:
    """
    Return genomic intervals for SNP awareness.

    - region_start: template origin for mapping hits (full VCF-spiked window when present).
    - query_start / query_end: 1-based inclusive Ensembl overlap interval, aligned with
      the same display window shown in the results UI (DISPLAY_FLANK).
    """
    genomic_pos = getattr(var_info, "genomic_pos", None)
    if not genomic_pos or not genomic_pos.get("pos"):
        return None

    chromosome = str(genomic_pos["chr"])
    pos = genomic_pos["pos"]
    if pos and isinstance(pos[0], list):
        starts = [int(p[0]) for p in pos]
        ends = [int(p[-1]) for p in pos]
        variant_start = min(starts)
        variant_end = max(ends)
    else:
        variant_start = int(pos[0])
        variant_end = int(pos[-1]) if len(pos) > 1 else int(pos[0])

    region_start_attr = getattr(var_info, "sequence_region_start", None)
    ref_seq = getattr(var_info, "ref_seq", None)
    if region_start_attr is not None and ref_seq:
        region_start = int(region_start_attr)
        region_end = region_start + len(ref_seq) - 1
    else:
        region_start = max(1, variant_start - flank)
        region_end = variant_end + flank

    rel = getattr(var_info, "relative_pos", None)
    display_start = 0
    display_end = (
        len(ref_seq) if ref_seq else (variant_end - variant_start + 1 + 2 * flank)
    )

    if rel and ref_seq:
        var_lo, var_hi = int(rel[0]), int(rel[1])
        if primer_target:
            target_start, target_len = primer_target
        else:
            target_start = max(0, var_lo - flank)
            target_len = (var_hi - var_lo + 1) + 2 * flank
        display_start, display_end = compute_display_bounds(
            len(ref_seq),
            var_lo,
            var_hi,
            target_start,
            target_len,
            primer_pairs or [],
        )
        query_start = region_start + display_start
        query_end = region_start + display_end - 1
    else:
        query_start = max(1, variant_start - flank)
        query_end = variant_end + flank

    return {
        "chromosome": chromosome,
        "region_start": region_start,
        "region_end": region_end,
        "query_start": max(1, query_start),
        "query_end": query_end,
        "display_start": display_start,
        "display_end": display_end,
        "variant_start": variant_start,
        "variant_end": variant_end,
    }


def _genomic_to_template_pos(genomic_1based: int, region_start: int) -> int:
    return genomic_1based - region_start


def _normalize_variation_hit(
    hit: dict, region_start: int, region_end: int
) -> Optional[dict]:
    try:
        g_start = int(hit["start"])
        g_end = int(hit["end"])
    except (KeyError, TypeError, ValueError):
        return None

    if g_end < region_start or g_start > region_end:
        return None

    template_start = _genomic_to_template_pos(g_start, region_start)
    template_end = _genomic_to_template_pos(g_end, region_start)
    alleles = hit.get("alleles") or []
    if isinstance(alleles, list):
        alleles_str = "/".join(str(a) for a in alleles if a is not None)
    else:
        alleles_str = str(alleles)

    return {
        "id": hit.get("id", ""),
        "source": hit.get("source", ""),
        "genomic_start": g_start,
        "genomic_end": g_end,
        "template_start": template_start,
        "template_end": template_end,
        "alleles": alleles_str,
        "consequence_type": hit.get("consequence_type", ""),
    }


def _is_user_variant_hit(hit: dict, var_info: AllelicVariantInfo) -> bool:
    """Exclude dbSNP entries overlapping the user's variant in template coordinates."""
    rel = getattr(var_info, "relative_pos", None)
    if not rel:
        return False
    return _intervals_overlap(
        hit["template_start"],
        hit["template_end"],
        int(rel[0]),
        int(rel[1]),
    )


def _classify_pair(
    pair: PrimerPairResult,
    hits: List[dict],
) -> Tuple[str, List[dict]]:
    left_start = pair.left_relPos_start
    left_end = pair.left_relPos_end
    right_start = pair.right_relPos_start
    right_end = pair.right_relPos_end
    if any(v is None for v in (left_start, left_end, right_start, right_end)):
        return SNP_STATUS_NONE, []

    conflicts: List[dict] = []
    caution_hits: List[dict] = []
    product_start = min(left_end, right_end) + 1
    product_end = max(left_start, right_start) - 1

    for hit in hits:
        ts, te = hit["template_start"], hit["template_end"]
        in_forward = _intervals_overlap(ts, te, left_start, left_end)
        in_reverse = _intervals_overlap(ts, te, right_start, right_end)
        if in_forward or in_reverse:
            entry = {**hit, "primer": "forward" if in_forward else "reverse"}
            if in_forward and in_reverse:
                entry["primer"] = "forward+reverse"
            conflicts.append(entry)
        elif product_start <= product_end and _intervals_overlap(
            ts, te, product_start, product_end
        ):
            caution_hits.append(hit)

    if conflicts:
        return SNP_STATUS_CONFLICT, conflicts
    if caution_hits:
        return SNP_STATUS_CAUTION, caution_hits
    return SNP_STATUS_NONE, []


def annotate_primer_pairs_with_snp_awareness(
    var_info: AllelicVariantInfo,
    primer_pairs: List[PrimerPairResult],
    reference_genome: str,
    *,
    enabled: bool = True,
    primer_target: Optional[Tuple[int, int]] = None,
) -> dict:
    """
    Query Ensembl for known variants in the design region and annotate each primer pair.

    Returns a JSON-serializable summary stored on DesignResultsSummary.snp_analysis_data.
    """
    base_summary: dict[str, Any] = {
        "enabled": enabled,
        "status": SNP_STATUS_SKIPPED,
        "reference_genome": reference_genome,
        "region": None,
        "variant_count": 0,
        "hits": [],
        "pairs_with_binding_conflict": 0,
        "maf_threshold": COMMON_VARIANT_MAF_THRESHOLD,
        "variant_set": GNOMAD_VARIANT_SET,
        "message": "",
    }

    if not enabled:
        base_summary["message"] = "SNP check was not requested."
        for pair in primer_pairs:
            pair.snp_status = SNP_STATUS_SKIPPED
            pair.snp_conflicts = []
        return base_summary

    if isinstance(var_info, SequenceVariantInfo) and not getattr(
        var_info, "genomic_pos", None
    ):
        base_summary["message"] = (
            "SNP check requires genomic or transcript coordinates "
            "(not available for sequence-only input)."
        )
        for pair in primer_pairs:
            pair.snp_status = SNP_STATUS_SKIPPED
            pair.snp_conflicts = []
        return base_summary

    region = get_design_region_genomic(
        var_info,
        primer_pairs=primer_pairs,
        primer_target=primer_target,
    )
    if not region:
        base_summary["message"] = "No genomic coordinates available for SNP check."
        for pair in primer_pairs:
            pair.snp_status = SNP_STATUS_SKIPPED
            pair.snp_conflicts = []
        return base_summary

    query_start = region["query_start"]
    query_end = region["query_end"]
    base_summary["region"] = {
        "chromosome": region["chromosome"],
        "start": query_start,
        "end": query_end,
    }

    try:
        client = EnsemblClient(ref_genome=reference_genome)
        raw_hits = client.get_overlapping_variations_for_region(
            region["chromosome"],
            query_start,
            query_end,
            variant_set=GNOMAD_VARIANT_SET,
        )
    except Exception as exc:
        LOGGER.exception("Ensembl variation overlap failed")
        base_summary["status"] = SNP_STATUS_ERROR
        base_summary["message"] = f"Ensembl SNP lookup failed: {exc}"
        for pair in primer_pairs:
            pair.snp_status = SNP_STATUS_ERROR
            pair.snp_conflicts = []
        return base_summary

    candidate_hits: List[dict] = []
    for raw in raw_hits:
        normalized = _normalize_variation_hit(
            raw, region["region_start"], region["region_end"]
        )
        if not normalized:
            continue
        if _is_user_variant_hit(normalized, var_info):
            continue
        candidate_hits.append(normalized)

    try:
        variation_details = client.get_variation_details_batch(
            [h["id"] for h in candidate_hits if h.get("id")]
        )
    except Exception as exc:
        LOGGER.exception("Ensembl variation frequency lookup failed")
        base_summary["status"] = SNP_STATUS_ERROR
        base_summary["message"] = f"Ensembl SNP frequency lookup failed: {exc}"
        for pair in primer_pairs:
            pair.snp_status = SNP_STATUS_ERROR
            pair.snp_conflicts = []
        return base_summary

    hits = _attach_maf_and_filter_common_variants(candidate_hits, variation_details)

    base_summary["hits"] = hits
    base_summary["variant_count"] = len(hits)
    base_summary["status"] = "ok"

    conflict_pairs = 0
    for pair in primer_pairs:
        status, details = _classify_pair(pair, hits)
        pair.snp_status = status
        pair.snp_conflicts = details
        if status == SNP_STATUS_CONFLICT:
            conflict_pairs += 1

    base_summary["pairs_with_binding_conflict"] = conflict_pairs
    pct = int(COMMON_VARIANT_MAF_THRESHOLD * 100)
    if hits:
        base_summary["message"] = (
            f"{len(hits)} common variant(s) in design region "
            f"(gnomAD via Ensembl, MAF > {pct}%). "
            f"{conflict_pairs} primer pair(s) overlap SNP binding sites."
        )
    else:
        base_summary["message"] = (
            f"No common variants in design region (gnomAD, MAF > {pct}%)."
        )

    return base_summary
