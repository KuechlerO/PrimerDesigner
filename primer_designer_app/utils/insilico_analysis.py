import os
import subprocess
from django.conf import settings
import logging
import tempfile
import json


LOGGER = logging.getLogger(__name__)

# Set all directory and file paths
REFERENCE_DIR = settings.REFERENCE_DATA_DIR
print(f"Reference directory set to: {REFERENCE_DIR}")
DICEY_DIR = os.path.join(settings.BASE_DIR, 'primer_designer_app', 'dicey_extras')
DICEY_CONFIG = os.path.join(DICEY_DIR, 'primer3_config')


def prepare_context_path(primer_settings) -> dict[str, str]:
    """Prepare the reference file path based on the requested genome."""
    rg = primer_settings.reference_genome
    dna_ref_file = f"Homo_sapiens.{rg}.dna.primary_assembly.fa.gz"

    if primer_settings.reference_genome == 'GRCh37':
        cdna_ref_file = 'gencode.v37lift37.transcripts.fa.gz'
    elif primer_settings.reference_genome == 'GRCh38':
        cdna_ref_file = 'gencode.v49.transcripts.fa.gz'
    else:
        raise ValueError(
            f"Unsupported reference genome: {primer_settings.reference_genome}"
        )

    LOGGER.debug(
        'Prepared reference paths: cdna=%s, dna=%s', cdna_ref_file, dna_ref_file
    )
    return {
        'cdna': os.path.join(REFERENCE_DIR, cdna_ref_file),
        'dna': os.path.join(REFERENCE_DIR, dna_ref_file),
    }


def insilico_reference_description(primer_settings) -> str:
    """
    Short sentence for the amplicon modal: which reference file and mode Dicey used.
    """
    if not getattr(primer_settings, 'do_insilico_pcr', False):
        return ''
    paths = prepare_context_path(primer_settings)
    assembly = primer_settings.reference_genome
    if primer_settings.context == 'genomic':
        fn = os.path.basename(paths['dna'])
        return (
            f"Genomic amplicons were identified using the reference genome file {fn} "
            f"({assembly})."
        )
    if primer_settings.context == 'transcriptomic':
        fn = os.path.basename(paths['cdna'])
        return (
            'Transcriptomic amplicons were identified using the transcriptome '
            f"reference file {fn} ({assembly})."
        )
    return ''


def run_dicey(primer_file, reference, output_gz, output_json):
    """Run Dicey and return the parsed JSON output."""
    dicey_cmd = [
        'dicey',
        'search',
        '-i',
        DICEY_CONFIG,
        '-l',
        '15000',
        '-o',
        output_gz,
        '-g',
        reference,
        primer_file,
    ]
    try:
        subprocess.run(dicey_cmd, check=True)
        subprocess.run(['gunzip', '-f', output_gz], check=True)
        with open(output_json, 'r', encoding='utf-8') as fh:
            return json.load(fh)
    except (subprocess.CalledProcessError, OSError, json.JSONDecodeError) as e:
        LOGGER.exception('Dicey failed or output could not be read: %s', e)
        return None


def process_primer_pair(pair, index, context_path, temp_dir) -> dict:
    """Process a single primer pair and return amplicons, seq, and status."""
    from primer_designer_app.utils.primer_utils import (
        INSILICO_ERROR,
        INSILICO_OK,
        INSILICO_OK_EMPTY,
    )

    LOGGER.debug(
        'Processing primer pair %d: Forward %s, Reverse %s',
        index,
        pair.left_seq,
        pair.right_seq,
    )
    primer_file = os.path.join(temp_dir, f"primer_{index}.fa")
    try:
        with open(primer_file, 'w', encoding='utf-8') as f:
            f.write(
                f">FGA_f{index}\n{pair.left_seq}\n>FGA_r{index}\n{pair.right_seq}\n"
            )
    except OSError as e:
        LOGGER.exception('Could not write primer file %s: %s', primer_file, e)
        detail = str(e)
        return {
            'amplicons': [],
            'insilico_seq': '',
            'insilico_status': INSILICO_ERROR,
            'insilico_error_detail': detail,
        }

    in_silico_output = run_dicey(
        primer_file,
        context_path,
        os.path.join(temp_dir, 'out.json.gz'),
        os.path.join(temp_dir, 'out.json'),
    )
    LOGGER.debug('Dicey output for primer pair %d: %s', index, in_silico_output)
    if not in_silico_output:
        return {
            'amplicons': [],
            'insilico_seq': '',
            'insilico_status': INSILICO_ERROR,
            'insilico_error_detail': 'Dicey failed or output could not be read',
        }

    amplicons = in_silico_output.get('data', {}).get('amplicons', [])
    insilico_seq = ''

    if len(amplicons) == 1:
        gen_pos_for = pair.left_relPos_start
        gen_pos_end = pair.right_relPos_end
        amp0 = amplicons[0]
        if gen_pos_for == amp0.get('ForPos') and gen_pos_end == amp0.get('RevEnd'):
            insilico_seq = amp0.get('Seq', '')

    if len(amplicons) == 0:
        status = INSILICO_OK_EMPTY
    else:
        status = INSILICO_OK

    return {
        'amplicons': amplicons,
        'insilico_seq': insilico_seq,
        'insilico_status': status,
        'insilico_error_detail': None,
    }


def do_insilico_analysis(primer_settings, primer_pairs: list) -> None:
    """
    Run Dicey in-silico search per primer pair and store results on
    each PrimerPairResult.
    """
    # Main logic
    reference_paths = prepare_context_path(primer_settings)
    if primer_settings.context == 'transcriptomic':
        context_path = reference_paths['cdna']
    elif primer_settings.context == 'genomic':
        context_path = reference_paths['dna']
    else:
        raise ValueError(
            f"Unsupported context for insilico-analysis: {primer_settings.context}"
        )

    with tempfile.TemporaryDirectory() as temp_dir:
        for i, pair in enumerate(primer_pairs):
            result = process_primer_pair(pair, i, context_path, temp_dir)
            pair.amplicons = result['amplicons']
            LOGGER.debug('Amplicons for pair %d: %s', i, pair.amplicons)
            pair.insilico_seq = result['insilico_seq']
            pair.insilico_status = result['insilico_status']
            pair.insilico_error_detail = result.get('insilico_error_detail')
