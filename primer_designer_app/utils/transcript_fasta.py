"""Local Ensembl cDNA/CDS FASTA access via samtools faidx."""

from __future__ import annotations

import logging
import os
import subprocess

from django.conf import settings

from primer_designer_app.exceptions import InvalidReferenceSequenceError

LOGGER = logging.getLogger(__name__)

_CDNA_FILES = {
    "GRCh37": "Homo_sapiens.GRCh37.cdna.all.fa.gz",
    "GRCh38": "Homo_sapiens.GRCh38.cdna.all.fa.gz",
}
_CDS_FILES = {
    "GRCh37": "Homo_sapiens.GRCh37.cds.all.fa.gz",
    "GRCh38": "Homo_sapiens.GRCh38.cds.all.fa.gz",
}


def transcript_fasta_path(ref_genome: str, seq_type: str) -> str:
    """Return absolute path to the local Ensembl cDNA or CDS FASTA for ``ref_genome``."""
    seq_type = (seq_type or "").lower()
    if seq_type == "cdna":
        mapping = _CDNA_FILES
    elif seq_type == "cds":
        mapping = _CDS_FILES
    else:
        raise ValueError(f"Unsupported seq_type for local transcript FASTA: {seq_type}")

    if ref_genome not in mapping:
        raise ValueError(
            f"Unsupported reference genome: {ref_genome}. "
            "Supported values are 'GRCh38' and 'GRCh37'."
        )

    reference_dir = settings.REFERENCE_DATA_DIR
    return os.path.join(reference_dir, mapping[ref_genome])


def _parse_faidx_stdout(stdout: str) -> str:
    """Strip FASTA header/newlines from ``samtools faidx`` stdout."""
    lines = []
    for line in (stdout or "").splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        lines.append(line)
    return "".join(lines)


def fetch_transcript_sequence(
    transcript_id: str, seq_type: str, ref_genome: str
) -> str:
    """Fetch one transcript sequence from the local faidx'd Ensembl FASTA.

    Requires an exact ID match (e.g. ``ENST00000424479.7``). Missing files or IDs
    raise ``InvalidReferenceSequenceError`` with setup-oriented guidance.
    """
    transcript_id = (transcript_id or "").strip()
    if not transcript_id:
        raise InvalidReferenceSequenceError("Transcript ID is required.")

    fasta_path = transcript_fasta_path(ref_genome, seq_type)
    fai_path = fasta_path + ".fai"

    if not os.path.isfile(fasta_path):
        raise InvalidReferenceSequenceError(
            f"Local {seq_type.upper()} reference file not found: {fasta_path}. "
            f"Download the Ensembl {ref_genome} {seq_type}.all FASTA into "
            "REFERENCE_DATA_DIR and run samtools faidx (see README)."
        )
    if not os.path.isfile(fai_path):
        raise InvalidReferenceSequenceError(
            f"Missing FASTA index for {fasta_path} (expected {fai_path}). "
            "Run: samtools faidx " + os.path.basename(fasta_path)
        )

    LOGGER.debug(
        "Fetching local %s sequence for %s from %s",
        seq_type,
        transcript_id,
        fasta_path,
    )
    try:
        completed = subprocess.run(
            ["samtools", "faidx", fasta_path, transcript_id],
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as exc:
        raise InvalidReferenceSequenceError(
            "samtools is not available on PATH; required to read local transcript FASTA."
        ) from exc

    if completed.returncode != 0:
        err = (completed.stderr or completed.stdout or "").strip()
        raise InvalidReferenceSequenceError(
            f"Transcript {transcript_id} not found in local {seq_type.upper()} FASTA "
            f"({os.path.basename(fasta_path)}). Ensure the FASTA matches {ref_genome} "
            f"and includes this Ensembl ID/version. samtools: {err or 'no details'}"
        )

    seq = _parse_faidx_stdout(completed.stdout)
    if not seq:
        raise InvalidReferenceSequenceError(
            f"Empty sequence returned for {transcript_id} from {fasta_path}."
        )
    return seq
