import os
import subprocess
from django.conf import settings
import logging
import tempfile
import json


LOGGER = logging.getLogger(__name__)

# Set all directory and file paths
REFERENCE_DIR = settings.REFERENCE_DATA_DIR
DICEY_DIR = os.path.join(settings.BASE_DIR, "primer_designer_app", "dicey_extras")
DICEY_CONFIG = os.path.join(DICEY_DIR, "primer3_config")


def prepare_context_path(primer_settings) -> dict[str, str]:
    """Prepare the reference file path based on the requested genome."""
    dna_ref_file = f"{primer_settings.reference_genome}.primary_assembly.genome.fa.gz"

    if primer_settings.reference_genome == "GRCh37":
        cdna_ref_file = "gencode.v37lift37.transcripts.fa.gz"
    elif primer_settings.reference_genome == "GRCh38":
        cdna_ref_file = "gencode.v49.transcripts.fa.gz"
    else:
        raise ValueError(f"Unsupported reference genome: {primer_settings.reference_genome}")
    
    return {
        "cdna": f"{REFERENCE_DIR}{cdna_ref_file}",
        "dna": f"{REFERENCE_DIR}{dna_ref_file}",
    }


def run_dicey(primer_file, reference, output_gz, output_json):
    """Run Dicey and return the parsed JSON output."""
    dicey_cmd = [
        "dicey",
        "search",
        "-i",
        DICEY_CONFIG,
        "-l",
        "15000",
        "-o",
        output_gz,
        "-g",
        reference,
        primer_file,
    ]
    try:
        subprocess.run(dicey_cmd, check=True)
        subprocess.run(["gunzip", "-f", output_gz], check=True)
        with open(output_json, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except (subprocess.CalledProcessError, OSError, json.JSONDecodeError) as e:
        LOGGER.exception("Dicey failed or output could not be read: %s", e)
        return None


def process_primer_pair(pair, index, context_path, temp_dir):
    """Process a single primer pair and return its results."""
    LOGGER.debug(
        "Processing primer pair %d: Forward %s, Reverse %s",
        index,
        pair.left_seq,
        pair.right_seq,
    )
    primer_file = os.path.join(temp_dir, f"primer_{index}.fa")
    try:
        with open(primer_file, "w", encoding="utf-8") as f:
            f.write(
                f">FGA_f{index}\n{pair.left_seq}\n>FGA_r{index}\n{pair.right_seq}\n"
            )
    except OSError as e:
        LOGGER.exception("Could not write primer file %s: %s", primer_file, e)
        return {"amplicons": [], "insilico_seq": ""}

    # Run Dicey for the primary reference
    in_silico_output = run_dicey(
        primer_file,
        context_path,
        os.path.join(temp_dir, "out.json.gz"),
        os.path.join(temp_dir, "out.json"),
    )
    LOGGER.debug("Dicey output for primer pair %d: %s", index, in_silico_output)
    if not in_silico_output:
        return {"amplicons": [], "insilico_seq": ""}

    amplicons = in_silico_output.get("data", {}).get("amplicons", [])
    insilico_seq = ""

    amplicons = amplicons
    if len(amplicons) == 1:
        gen_pos_for = pair.left_relPos_start
        gen_pos_end = pair.right_relPos_end
        amp0 = amplicons[0]
        if gen_pos_for == amp0.get("ForPos") and gen_pos_end == amp0.get("RevEnd"):
            insilico_seq = amp0.get("Seq", "")

    return {
        "amplicons": amplicons,
        "insilico_seq": insilico_seq,
    }


def do_insilico_analysis(primer_settings, primer_pairs: list) -> None:
    """
    Run Dicey in-silico search per primer pair and store results on
    each PrimerPairResult.
    """
    # Main logic
    reference_paths = prepare_context_path(primer_settings)
    if primer_settings.context == "transcriptomic":
        context_path = reference_paths["cdna"]
    elif primer_settings.context == "genomic":
        context_path = reference_paths["dna"]
    else:
        LOGGER.error("Invalid context '%s' in primer settings", primer_settings.context)
        return

    with tempfile.TemporaryDirectory() as temp_dir:
        for i, pair in enumerate(primer_pairs):
            result = process_primer_pair(pair, i, context_path, temp_dir)
            pair.amplicons = result["amplicons"]
            pair.insilico_seq = result["insilico_seq"]
