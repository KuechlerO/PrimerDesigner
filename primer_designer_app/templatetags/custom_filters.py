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

register = template.Library()


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
    return extract_amplicon_summary(amplicon_dict)


@register.filter
def amplicon_chrom_display(amplicon_dict):
    return amplicon_chrom_label(amplicon_dict)


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
    return format_penalty_score(value)
