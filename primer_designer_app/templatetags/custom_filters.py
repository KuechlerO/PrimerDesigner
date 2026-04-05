from django import template

register = template.Library()


def _chrom_parts(amplicon_dict):
    chrom = amplicon_dict.get("Chrom") or ""
    return chrom.split("|")


@register.filter
def extract_amplicon_info(amplicon_dict, delimiter='|'):
    parts = _chrom_parts(amplicon_dict)
    if len(parts) == 1:
        return f"{amplicon_dict['Chrom']}: {amplicon_dict['ForPos']} - {amplicon_dict['RevEnd']}"
    elif len(parts) > 1:
        transcript_id = parts[0]
        gene_symbol = parts[-4]

        return f"{gene_symbol} ({transcript_id}): {amplicon_dict['ForPos']} - {amplicon_dict['RevEnd']}"
    else:
        return "Invalid amplicon information"


@register.filter
def amplicon_chrom_display(amplicon_dict):
    """Genomic: chromosome name. Transcriptome: gene (ENST…) only, matching summary semantics."""
    parts = _chrom_parts(amplicon_dict)
    if len(parts) == 1:
        return parts[0]
    if len(parts) > 1:
        transcript_id = parts[0]
        gene_symbol = parts[-4]
        return f"{gene_symbol} ({transcript_id})"
    return ""


@register.filter
def penalty_two_decimals(value):
    if value is None or value == "":
        return ""
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return value