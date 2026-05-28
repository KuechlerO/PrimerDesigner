"""Parse VCF files and spike variants into a reference sequence window."""

from __future__ import annotations

import gzip
import io
import logging
from dataclasses import dataclass
from typing import BinaryIO, List, Optional, Tuple

LOGGER = logging.getLogger(__name__)

MAX_VCF_BYTES = 10 * 1024 * 1024
MAX_VCF_RECORDS = 5000


@dataclass
class VcfRecord:
    """One VCF row (first ALT allele used when spiking)."""

    chrom: str
    pos: int  # 1-based
    ref: str
    alt: str
    rsid: str = ""

    @property
    def end(self) -> int:
        return self.pos + len(self.ref) - 1


@dataclass
class AppliedVcfVariant:
    """Metadata for variants successfully spiked into the template."""

    id: str
    chrom: str
    pos: int
    ref: str
    alt: str
    template_start: int
    template_end: int


def normalize_chromosome(chrom: str) -> str:
    """Normalize to Ensembl-style chromosome (no chr prefix, MT→M)."""
    c = (chrom or "").strip()
    if not c:
        return ""
    lower = c.lower()
    if lower.startswith("chr"):
        c = c[3:]
    if c == "MT":
        c = "M"
    return c.upper()


def _open_vcf_stream(uploaded_file) -> BinaryIO:
    name = (getattr(uploaded_file, "name", "") or "").lower()
    raw = uploaded_file.read(MAX_VCF_BYTES + 1)
    if len(raw) > MAX_VCF_BYTES:
        raise ValueError(
            f"VCF file exceeds maximum size ({MAX_VCF_BYTES // (1024 * 1024)} MB)."
        )
    if name.endswith(".gz"):
        return gzip.open(io.BytesIO(raw), "rt", encoding="utf-8", errors="replace")
    return io.StringIO(raw.decode("utf-8", errors="replace"))


def parse_vcf_upload(uploaded_file, target_chrom: str) -> List[VcfRecord]:
    """
    Parse an uploaded VCF/VCF.GZ; return records on target_chrom (normalized).
    Skips structural alleles and non-ACGTN ref/alt.
    """
    if not uploaded_file:
        return []

    target = normalize_chromosome(target_chrom)
    records: List[VcfRecord] = []
    stream = _open_vcf_stream(uploaded_file)

    try:
        for line in stream:
            if isinstance(line, bytes):
                line = line.decode("utf-8", errors="replace")
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            chrom = normalize_chromosome(parts[0])
            if chrom != target:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            rsid = parts[2] if parts[2] != "." else ""
            ref = parts[3].upper()
            alts = parts[4].split(",")
            if not ref or not alts or alts[0] == ".":
                continue
            alt = alts[0].upper()
            if alt == "*" or alt.startswith("<") or ref.startswith("<"):
                continue
            if not _is_nucleotide_allele(ref) or not _is_nucleotide_allele(alt):
                continue
            records.append(VcfRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, rsid=rsid))
            if len(records) >= MAX_VCF_RECORDS:
                LOGGER.warning("VCF truncated at %s records", MAX_VCF_RECORDS)
                break
    finally:
        if hasattr(stream, "close"):
            stream.close()

    records.sort(key=lambda r: (r.pos, len(r.ref)))
    return records


def _is_nucleotide_allele(allele: str) -> bool:
    return bool(allele) and set(allele) <= set("ACGTN")


def compute_fetch_window(
    primary_start: int,
    primary_end: int,
    vcf_records: List[VcfRecord],
    flank: int,
) -> Tuple[int, int]:
    """1-based inclusive region to fetch from Ensembl."""
    starts = [primary_start]
    ends = [primary_end]
    for rec in vcf_records:
        starts.append(rec.pos)
        ends.append(rec.end)
    region_start = max(1, min(starts) - flank)
    region_end = max(ends) + flank
    return region_start, region_end


def _offset_before(genomic_pos: int, applied: List[Tuple[int, int, int]]) -> int:
    """Sum of (len(alt)-len(ref)) for spikes with POS < genomic_pos."""
    delta = 0
    for pos, ref_len, alt_len in applied:
        if pos < genomic_pos:
            delta += alt_len - ref_len
    return delta


def template_range_for_genomic(
    region_start_1based: int,
    genomic_start: int,
    genomic_end: int,
    applied: List[Tuple[int, int, int]],
) -> Tuple[int, int]:
    """0-based inclusive template indices after VCF spikes (before user primary edit)."""
    offset_start = _offset_before(genomic_start, applied)
    offset_end = _offset_before(genomic_end, applied)
    t_start = genomic_start - region_start_1based + offset_start
    t_end = genomic_end - region_start_1based + offset_end
    return max(0, t_start), max(0, t_end)


def spike_vcf_variants(
    ref_seq: str,
    region_start_1based: int,
    records: List[VcfRecord],
    *,
    skip_interval: Optional[Tuple[int, int]] = None,
) -> Tuple[str, List[AppliedVcfVariant], List[Tuple[int, int, int]]]:
    """
    Apply VCF records in ascending genomic order onto ref_seq.

    skip_interval: optional 1-based inclusive primary interval — overlapping VCF rows
    are skipped so the user's primary variant is applied separately.

    Returns (mutated_sequence, applied_metadata, spike_deltas for coord mapping).
    """
    seq = ref_seq
    applied_meta: List[AppliedVcfVariant] = []
    deltas: List[Tuple[int, int, int]] = []

    for rec in records:
        if skip_interval:
            p_start, p_end = skip_interval
            if rec.pos <= p_end and rec.end >= p_start:
                LOGGER.debug(
                    "Skipping VCF %s at %s:%s (overlaps primary variant)",
                    rec.rsid or ".",
                    rec.chrom,
                    rec.pos,
                )
                continue

        offset = _offset_before(rec.pos, deltas)
        idx = rec.pos - region_start_1based + offset
        if idx < 0 or idx + len(rec.ref) > len(seq):
            LOGGER.warning(
                "Skipping VCF %s:%s — outside fetched window",
                rec.chrom,
                rec.pos,
            )
            continue
        ref_in_seq = seq[idx : idx + len(rec.ref)]
        if ref_in_seq.upper() != rec.ref.upper():
            LOGGER.warning(
                "Skipping VCF %s:%s — REF mismatch (VCF %s, reference %s)",
                rec.chrom,
                rec.pos,
                rec.ref,
                ref_in_seq,
            )
            continue

        seq = seq[:idx] + rec.alt + seq[idx + len(rec.ref) :]
        t_end = idx + len(rec.alt) - 1
        applied_meta.append(
            AppliedVcfVariant(
                id=rec.rsid or f"{rec.chrom}:{rec.pos}",
                chrom=rec.chrom,
                pos=rec.pos,
                ref=rec.ref,
                alt=rec.alt,
                template_start=idx,
                template_end=t_end,
            )
        )
        deltas.append((rec.pos, len(rec.ref), len(rec.alt)))

    return seq, applied_meta, deltas
