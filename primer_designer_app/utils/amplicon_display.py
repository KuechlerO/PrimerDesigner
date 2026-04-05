"""Shared formatting for in-silico amplicons (UI templates and Word reports)."""


def _chrom_parts(amplicon_dict: dict) -> list[str]:
    chrom = amplicon_dict.get('Chrom') or ''
    return chrom.split('|')


def extract_amplicon_summary(amplicon_dict: dict) -> str:
    """Same semantics as template filter extract_amplicon_info."""
    parts = _chrom_parts(amplicon_dict)
    if len(parts) == 1:
        return (
            f"{amplicon_dict['Chrom']}: {amplicon_dict['ForPos']} - "
            f"{amplicon_dict['RevEnd']}"
        )
    if len(parts) > 1:
        transcript_id = parts[0]
        gene_symbol = parts[-4]
        return (
            f"{gene_symbol} ({transcript_id}): {amplicon_dict['ForPos']} - "
            f"{amplicon_dict['RevEnd']}"
        )
    return 'Invalid amplicon information'


def amplicon_chrom_label(amplicon_dict: dict) -> str:
    """Genomic: chromosome. Transcriptome: gene (ENST)."""
    parts = _chrom_parts(amplicon_dict)
    if len(parts) == 1:
        return parts[0]
    if len(parts) > 1:
        transcript_id = parts[0]
        gene_symbol = parts[-4]
        return f"{gene_symbol} ({transcript_id})"
    return ''


def format_penalty_score(value) -> str:
    if value is None or value == '':
        return ''
    try:
        return f"{float(value):.2f}"
    except (TypeError, ValueError):
        return str(value)


def truncate_product_seq(seq, max_len: int = 64) -> str:
    if seq is None:
        return ''
    s = str(seq)
    if len(s) <= max_len:
        return s
    return s[: max_len - 1] + '…'
