from django import template

from primer_designer_app.utils.amplicon_display import (
    amplicon_chrom_label,
    extract_amplicon_summary,
    format_penalty_score,
)
from primer_designer_app.utils.primer_utils import (
    INSILICO_ERROR,
    INSILICO_NOT_APPLICABLE,
    INSILICO_OK,
    INSILICO_OK_EMPTY,
)
from primer_designer_app.utils.snp_awareness import (
    SNP_STATUS_CAUTION,
    SNP_STATUS_CONFLICT,
    SNP_STATUS_ERROR,
    SNP_STATUS_NONE,
    SNP_STATUS_SKIPPED,
)

register = template.Library()


@register.filter
def insilico_cell_class(status):
    """CSS classes for the in-silico / amplicons table cell."""
    return {
        INSILICO_OK: "amplicon-cell amplicon-cell--ok",
        INSILICO_OK_EMPTY: "amplicon-cell amplicon-cell--empty",
        INSILICO_NOT_APPLICABLE: "amplicon-cell amplicon-cell--na",
        INSILICO_ERROR: "amplicon-cell amplicon-cell--error",
    }.get(status or "", "amplicon-cell amplicon-cell--unknown")


@register.filter
def extract_amplicon_info(amplicon_dict, delimiter="|"):
    return extract_amplicon_summary(amplicon_dict)


@register.filter
def amplicon_chrom_display(amplicon_dict):
    return amplicon_chrom_label(amplicon_dict)


@register.filter
def insilico_ok_variant_class(pair):
    """When status is ok: highlight single vs multiple amplicons."""
    if getattr(pair, "insilico_status", None) != INSILICO_OK:
        return ""
    n = len(pair.amplicons or [])
    if n == 1:
        return "amplicon-cell--ok-single"
    return "amplicon-cell--ok-multi"


@register.filter
def penalty_two_decimals(value):
    return format_penalty_score(value)


@register.filter
def snp_cell_class(status):
    return {
        SNP_STATUS_NONE: "snp-cell snp-cell--none",
        SNP_STATUS_CAUTION: "snp-cell snp-cell--caution",
        SNP_STATUS_CONFLICT: "snp-cell snp-cell--conflict",
        SNP_STATUS_SKIPPED: "snp-cell snp-cell--skipped",
        SNP_STATUS_ERROR: "snp-cell snp-cell--error",
    }.get(status or "", "snp-cell snp-cell--na")


@register.filter
def snp_status_label(status):
    return {
        SNP_STATUS_NONE: "Clear",
        SNP_STATUS_CAUTION: "Caution",
        SNP_STATUS_CONFLICT: "Conflict",
        SNP_STATUS_SKIPPED: "N/A",
        SNP_STATUS_ERROR: "Error",
    }.get(status or "", "—")


@register.filter
def snp_conflict_summary(conflicts):
    if not conflicts:
        return ""
    parts = []
    for hit in conflicts[:3]:
        rsid = hit.get("id", "")
        primer = hit.get("primer", "")
        alleles = hit.get("alleles", "")
        parts.append(f"{rsid} ({primer}, {alleles})")
    text = "; ".join(parts)
    if len(conflicts) > 3:
        text += f" (+{len(conflicts) - 3} more)"
    return text
