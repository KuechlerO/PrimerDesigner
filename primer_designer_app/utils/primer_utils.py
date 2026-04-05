import logging
from dataclasses import dataclass, asdict
from typing import List, Optional

import primer3
from primer_designer_app.utils.variant_info import (
    VariantInfo,
    TranscriptVariantInfo,
    GenomicVariantInfo,
    VARIANT_FLANKING,
)
from primer_designer_app.utils.insilico_analysis import do_insilico_analysis

LOGGER = logging.getLogger(__name__)

# In-silico (Dicey) outcome per primer pair — persisted on PrimerPairResult
INSILICO_NOT_APPLICABLE = 'not_applicable'
INSILICO_ERROR = 'error'
INSILICO_OK_EMPTY = 'ok_empty'
INSILICO_OK = 'ok'


@dataclass
class PrimerPairResult:
    # Container for dicey primer pair result
    index: int
    left_seq: str
    right_seq: str
    penalty: float
    product_size: int
    product_tm: Optional[float] = None
    left_relPos_start: Optional[int] = None  # 0-based position
    left_relPos_end: Optional[int] = None  # 0-based position
    right_relPos_start: Optional[int] = None  # 0-based position
    right_relPos_end: Optional[int] = None  # 0-based position
    gc: Optional[List[float]] = None
    tm: Optional[List[float]] = None
    insilico_seq: Optional[str] = ''
    amplicons: Optional[List[str]] = None
    insilico_status: Optional[str] = None
    insilico_error_detail: Optional[str] = None

    def to_dict(self):
        return asdict(self)


class PrimerSearchResults:
    def __init__(self, primer3_obj=None):
        # grouped per-pair results
        self.primer_pairs: List[PrimerPairResult] = []
        if primer3_obj:
            # primer3_obj is the dict returned by primer3.bindings.design_primers
            self.primer_pairs = get_primers_from_primer3(primer3_obj)

        self.mapped_primer_positions = {'primerF_starts': [], 'primerR_ends': []}

    @classmethod
    def from_dict(cls, data: dict):
        """Create a PrimerSearchResults object from a dictionary."""
        instance = cls()
        instance.mapped_primer_positions = data.get(
            'mapped_primer_positions', {'primerF_starts': [], 'primerR_ends': []}
        )

        instance.primer_pairs = [
            primer_pair_from_dict(pair) for pair in data.get('primer_pairs', [])
        ]
        return instance

    def to_dict(self):
        return asdict(self)

    def load_primer_start_and_end_pos(self, varInfo: VariantInfo):
        """Calculate genomic positions of primers based on variant info."""
        gene_pos = varInfo.get_genomic_pos()[0]
        # clear existing mapped positions before filling
        self.mapped_primer_positions['primerF_starts'].clear()
        self.mapped_primer_positions['primerR_ends'].clear()

        for pair in self.primer_pairs:
            if pair.left_relPos_start is None or pair.right_relPos_end is None:
                # skip malformed pair
                self.mapped_primer_positions['primerF_starts'].append(None)
                self.mapped_primer_positions['primerR_ends'].append(None)
                continue
            self.mapped_primer_positions['primerF_starts'].append(
                gene_pos - VARIANT_FLANKING + pair.left_relPos_start
            )
            self.mapped_primer_positions['primerR_ends'].append(
                gene_pos - VARIANT_FLANKING + pair.right_relPos_end
            )
        return self.mapped_primer_positions


def _infer_legacy_insilico_status(pair_dict: dict) -> str:
    """Best-effort status for JSON saved before insilico_status existed."""
    amps = pair_dict.get('amplicons')
    if amps is None:
        return INSILICO_OK_EMPTY
    if amps == ['N/A'] or (
        isinstance(amps, list) and len(amps) == 1 and amps[0] == 'N/A'
    ):
        return INSILICO_NOT_APPLICABLE
    if isinstance(amps, list) and len(amps) > 0 and isinstance(amps[0], dict):
        return INSILICO_OK
    if amps == []:
        return INSILICO_OK_EMPTY
    return INSILICO_OK_EMPTY


def primer_pair_from_dict(pair: dict) -> PrimerPairResult:
    d = dict(pair)
    if d.get('insilico_status') is None:
        d['insilico_status'] = _infer_legacy_insilico_status(d)
    d.setdefault('insilico_error_detail', None)
    return PrimerPairResult(**d)


def get_primers_from_primer3(dicey_primer) -> List['PrimerPairResult']:
    """
    Create a list of PrimerPairResult from a primer3/dicey result dict.
    :param dicey_primer: dict returned by primer3.bindings.design_primers
    :return: List of PrimerPairResult (ordered by penalty)
    """
    pairs: List[PrimerPairResult] = []
    pair_count = int(dicey_primer.get('PRIMER_PAIR_NUM_RETURNED', 0))
    for i in range(pair_count):
        left_key = f"PRIMER_LEFT_{i}"
        right_key = f"PRIMER_RIGHT_{i}"
        left_seq_key = f"{left_key}_SEQUENCE"
        right_seq_key = f"{right_key}_SEQUENCE"
        left_gc_key = f"{left_key}_GC_PERCENT"
        right_gc_key = f"{right_key}_GC_PERCENT"
        left_tm_key = f"{left_key}_TM"
        right_tm_key = f"{right_key}_TM"
        pair_score_key = f"PRIMER_PAIR_{i}_PENALTY"
        pair_product_size_key = f"PRIMER_PAIR_{i}_PRODUCT_SIZE"
        product_tm_key = f"PRIMER_PAIR_{i}_PRODUCT_TM"

        if left_seq_key in dicey_primer and right_seq_key in dicey_primer:
            try:
                left_pos = dicey_primer[left_key]
                right_pos = dicey_primer[right_key]
                left_relPos_start = int(left_pos[0])
                left_len = int(left_pos[1])
                right_relPos_start = int(right_pos[0])
                right_len = int(right_pos[1])
            except Exception:
                # skip malformed entries
                continue

            pair = PrimerPairResult(
                index=i,
                left_seq=dicey_primer[left_seq_key],
                right_seq=dicey_primer[right_seq_key],
                penalty=round(float(dicey_primer.get(pair_score_key, 0.0)), 2),
                product_size=int(dicey_primer.get(pair_product_size_key, 0)),
                product_tm=round(float(dicey_primer.get(product_tm_key, 0.0)), 1),
                # Position is 0-based!
                left_relPos_start=left_relPos_start,
                left_relPos_end=left_relPos_start + (left_len - 1),
                right_relPos_start=right_relPos_start - (right_len - 1),
                right_relPos_end=right_relPos_start,
                gc=[
                    round(float(dicey_primer.get(left_gc_key, 0.0)), 1),
                    round(float(dicey_primer.get(right_gc_key, 0.0)), 1),
                ],
                tm=[
                    round(float(dicey_primer.get(left_tm_key, 0.0)), 1),
                    round(float(dicey_primer.get(right_tm_key, 0.0)), 1),
                ],
                amplicons=[],
                insilico_seq='',
                insilico_status=None,
                insilico_error_detail=None,
            )
            pairs.append(pair)

    # Ordered by penalty (lowest first)
    return pairs


def primer3_design_primers(
    primSet_obj, varInfo_obj: VariantInfo
) -> PrimerSearchResults:
    opt_tm = primSet_obj.tm
    LOGGER.info(f"Designing primers with target: {primSet_obj.target}")
    LOGGER.debug(f"Variant info object sequence: {varInfo_obj.get_seq("mutated")}")

    # Call primer3 to design primers
    primer3_obj = primer3.bindings.design_primers(
        seq_args={
            'SEQUENCE_ID': 'dummy_id',  # Sequence ID
            'SEQUENCE_TEMPLATE': varInfo_obj.get_seq(
                'mutated'
            ),  # full mutated sequence
            'SEQUENCE_TARGET': primSet_obj.target,  # [80,100]
        },
        # TODO: All these settings should be adjustable by the user
        global_args={
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,  # internal oligo
            'PRIMER_INTERNAL_MAX_SELF_END': 8,  # internal oligo -> TODO: Remove?!
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 22,
            'PRIMER_OPT_TM': opt_tm,
            'PRIMER_MIN_TM': opt_tm - 2,
            'PRIMER_MAX_TM': opt_tm + 2,
            'PRIMER_OPT_GC_PERCENT': primSet_obj.gc,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_GC_CLAMP': 1,
            'PRIMER_MAX_POLY_X': primSet_obj.max_poly_x,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,  # internal oligo -> TODO: Remove?! (default: 5)
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': primSet_obj.productsize_range,  # [100,500]
            'PRIMER_INSIDE_PENALTY': 1.0,  # dont allow primers inside the target (default 0 -> favors primers overlapping the target)
            # 'PRIMER_OUTSIDE_PENALTY': 0.0,
        },
    )

    LOGGER.debug(f"Primer3 output: {primer3_obj}")

    # Extract primer information
    prim3_res = PrimerSearchResults(primer3_obj=primer3_obj)

    LOGGER.debug(f"Using context: {primSet_obj.context} for in-silico analysis")
    # run in-silico analysis and populate amplicon summary
    if isinstance(varInfo_obj, (TranscriptVariantInfo, GenomicVariantInfo)):
        LOGGER.info('Loading primer start and end positions for genomic context')
        prim3_res.load_primer_start_and_end_pos(varInfo_obj)
        # run in-silico per primer pair
        LOGGER.info('Running in-silico analysis for designed primer pairs')
        do_insilico_analysis(primSet_obj, prim3_res.primer_pairs)
        LOGGER.info('In-silico analysis completed: ', prim3_res.primer_pairs)
    else:
        LOGGER.debug(
            'No genomic/transcript variant info; in-silico search not applicable'
        )
        for pair in prim3_res.primer_pairs:
            pair.amplicons = []
            pair.insilico_status = INSILICO_NOT_APPLICABLE
            pair.insilico_error_detail = None

    LOGGER.info(f"Primer positions: {prim3_res.mapped_primer_positions}")
    LOGGER.info(
        f"Primer positions relative to variant: {[(res.right_relPos_start, res.right_relPos_end) for res in prim3_res.primer_pairs]}"
    )
    return prim3_res
