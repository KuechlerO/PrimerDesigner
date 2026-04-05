from django import template

from primer_designer_app.utils.primer_utils import (
    INSILICO_ERROR,
    INSILICO_NOT_APPLICABLE,
    INSILICO_OK,
    INSILICO_OK_EMPTY,
)

register = template.Library()


def _chrom_parts(amplicon_dict):
    chrom = amplicon_dict.get('Chrom') or ''
    return chrom.split('|')


@register.filter
def insilico_cell_class(status):
    """CSS classes for the in-silico / amplicons table cell."""
    return {
        INSILICO_OK: 'amplicon-cell amplicon-cell--ok',
        INSILICO_OK_EMPTY: 'amplicon-cell amplicon-cell--empty',
        INSILICO_NOT_APPLICABLE: 'amplicon-cell amplicon-cell--na',
        INSILICO_ERROR: 'amplicon-cell amplicon-cell--error',
    }.get(status or '', 'amplicon-cell amplicon-cell--unknown')


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
        return 'Invalid amplicon information'


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
    return ''


@register.filter
def insilico_ok_variant_class(pair):
    """When status is ok: highlight single vs multiple amplicons."""
    if getattr(pair, 'insilico_status', None) != INSILICO_OK:
        return ''
    n = len(pair.amplicons or [])
    if n == 1:
        return 'amplicon-cell--ok-single'
    return 'amplicon-cell--ok-multi'


@register.filter
def penalty_two_decimals(value):
    if value is None or value == '':
        return ''
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return value
