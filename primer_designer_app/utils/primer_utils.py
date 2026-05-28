import logging
from dataclasses import dataclass, asdict
from typing import List, Optional

try:
    import primer3  # type: ignore
except ModuleNotFoundError:  # pragma: no cover
    primer3 = None
from primer_designer_app.utils.variant_info import (
    AllelicVariantInfo,
    IndelType,
    PrimerDesignSequence,
    StructuralVariantWindow,
    TranscriptVariantInfo,
    GenomicVariantInfo,
    SequenceVariantInfo,
    VARIANT_FLANKING,
)
from primer_designer_app.utils.insilico_analysis import do_insilico_analysis

LOGGER = logging.getLogger(__name__)

# In-silico (Dicey) outcome per primer pair — persisted on PrimerPairResult
INSILICO_NOT_APPLICABLE = "not_applicable"
INSILICO_ERROR = "error"
INSILICO_OK_EMPTY = "ok_empty"
INSILICO_OK = "ok"


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
    insilico_seq: Optional[str] = ""
    amplicons: Optional[List[str]] = None
    insilico_status: Optional[str] = None
    insilico_error_detail: Optional[str] = None
    snp_status: Optional[str] = None
    snp_conflicts: Optional[List[dict]] = None

    def to_dict(self):
        return asdict(self)


class PrimerSearchResults:
    def __init__(self, primer3_obj=None):
        # grouped per-pair results
        self.primer_pairs: List[PrimerPairResult] = []
        self.primer3_obj = primer3_obj
        if primer3_obj:
            # primer3_obj is the dict returned by primer3.bindings.design_primers
            self.primer_pairs = get_primers_from_primer3(primer3_obj)

        self.mapped_primer_positions = {"primerF_starts": [], "primerR_ends": []}

    @classmethod
    def from_dict(cls, data: dict):
        """Create a PrimerSearchResults object from a dictionary."""
        instance = cls()
        instance.mapped_primer_positions = data.get(
            "mapped_primer_positions", {"primerF_starts": [], "primerR_ends": []}
        )

        instance.primer_pairs = [
            primer_pair_from_dict(pair) for pair in data.get("primer_pairs", [])
        ]
        return instance

    def to_dict(self):
        # PrimerSearchResults is not a dataclass; serialize explicitly
        d = dict(self.__dict__)
        # The raw Primer3 result dict may contain bytes depending on environment/build;
        # it is not required for persistence and breaks JSON serialization.
        d.pop("primer3_obj", None)
        if "primer_pairs" in d:
            d["primer_pairs"] = [p.to_dict() for p in d.get("primer_pairs") or []]
        return d

    def load_primer_start_and_end_pos(self, varInfo: AllelicVariantInfo):
        """Calculate genomic positions of primers based on variant info."""
        gene_pos = varInfo.get_genomic_pos()[0]
        # clear existing mapped positions before filling
        self.mapped_primer_positions["primerF_starts"].clear()
        self.mapped_primer_positions["primerR_ends"].clear()

        for pair in self.primer_pairs:
            if pair.left_relPos_start is None or pair.right_relPos_end is None:
                # skip malformed pair
                self.mapped_primer_positions["primerF_starts"].append(None)
                self.mapped_primer_positions["primerR_ends"].append(None)
                continue
            self.mapped_primer_positions["primerF_starts"].append(
                gene_pos - VARIANT_FLANKING + pair.left_relPos_start
            )
            self.mapped_primer_positions["primerR_ends"].append(
                gene_pos - VARIANT_FLANKING + pair.right_relPos_end
            )
        return self.mapped_primer_positions


def _infer_legacy_insilico_status(pair_dict: dict) -> str:
    """Best-effort status for JSON saved before insilico_status existed."""
    amps = pair_dict.get("amplicons")
    if amps is None:
        return INSILICO_OK_EMPTY
    if amps == ["N/A"] or (
        isinstance(amps, list) and len(amps) == 1 and amps[0] == "N/A"
    ):
        return INSILICO_NOT_APPLICABLE
    if isinstance(amps, list) and len(amps) > 0 and isinstance(amps[0], dict):
        return INSILICO_OK
    if amps == []:
        return INSILICO_OK_EMPTY
    return INSILICO_OK_EMPTY


def primer_pair_from_dict(pair: dict) -> PrimerPairResult:
    d = dict(pair)
    if d.get("insilico_status") is None:
        d["insilico_status"] = _infer_legacy_insilico_status(d)
    d.setdefault("insilico_error_detail", None)
    d.setdefault("snp_status", None)
    d.setdefault("snp_conflicts", None)
    return PrimerPairResult(**d)


def get_primers_from_primer3(dicey_primer) -> List["PrimerPairResult"]:
    """
    Create a list of PrimerPairResult from a primer3/dicey result dict.
    :param dicey_primer: dict returned by primer3.bindings.design_primers
    :return: List of PrimerPairResult (ordered by penalty)
    """
    pairs: List[PrimerPairResult] = []
    pair_count = int(dicey_primer.get("PRIMER_PAIR_NUM_RETURNED", 0))
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
                right_pos_5p = int(right_pos[0])
                right_len = int(right_pos[1])
                right_relPos_start = right_pos_5p - (right_len - 1)
                right_relPos_end = right_pos_5p
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
                right_relPos_start=right_relPos_start,
                right_relPos_end=right_relPos_end,
                gc=[
                    round(float(dicey_primer.get(left_gc_key, 0.0)), 1),
                    round(float(dicey_primer.get(right_gc_key, 0.0)), 1),
                ],
                tm=[
                    round(float(dicey_primer.get(left_tm_key, 0.0)), 1),
                    round(float(dicey_primer.get(right_tm_key, 0.0)), 1),
                ],
                amplicons=[],
                insilico_seq="",
                insilico_status=None,
                insilico_error_detail=None,
            )
            pairs.append(pair)

    # Ordered by penalty (lowest first)
    return pairs


def build_primer3_global_args(prim_set) -> dict:
    """
    Build Primer3 global_args from model fields, then apply primer3_overrides (advanced POST keys).
    """
    opt_tm = float(prim_set.tm)
    base = {
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_INTERNAL_MAX_SELF_END": 8,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 22,
        "PRIMER_OPT_TM": opt_tm,
        "PRIMER_MIN_TM": opt_tm - 2,
        "PRIMER_MAX_TM": opt_tm + 2,
        "PRIMER_OPT_GC_PERCENT": float(prim_set.gc),
        "PRIMER_MIN_GC": 20.0,
        "PRIMER_MAX_GC": 80.0,
        "PRIMER_GC_CLAMP": 1,
        "PRIMER_MAX_POLY_X": prim_set.max_poly_x,
        "PRIMER_INTERNAL_MAX_POLY_X": 100,
        "PRIMER_SALT_MONOVALENT": 50.0,
        "PRIMER_DNA_CONC": 50.0,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_SELF_ANY": 12,
        "PRIMER_MAX_SELF_END": 8,
        "PRIMER_PAIR_MAX_COMPL_ANY": 12,
        "PRIMER_PAIR_MAX_COMPL_END": 8,
        "PRIMER_PRODUCT_SIZE_RANGE": prim_set.productsize_range,
        "PRIMER_INSIDE_PENALTY": 1.0,
    }
    overrides = getattr(prim_set, "primer3_overrides", None) or {}
    merged = {**base, **overrides}
    return merged


def primer3_design_primers(
    primSet_obj, varInfo_obj: PrimerDesignSequence
) -> PrimerSearchResults:
    if primer3 is None:
        raise RuntimeError("primer3 is not installed in this environment.")
    use_sequence_target = not isinstance(varInfo_obj, StructuralVariantWindow)
    LOGGER.info(
        "Designing primers (SEQUENCE_TARGET=%s, target=%s)",
        use_sequence_target,
        primSet_obj.target if use_sequence_target else None,
    )
    LOGGER.debug(
        "Variant info object sequence: %s",
        varInfo_obj.get_seq("mutated"),
    )

    global_args = build_primer3_global_args(primSet_obj)

    seq_args = {
        "SEQUENCE_ID": "dummy_id",
        "SEQUENCE_TEMPLATE": varInfo_obj.get_seq("mutated"),
    }
    if use_sequence_target:
        seq_args["SEQUENCE_TARGET"] = primSet_obj.target

    primer3_obj = primer3.bindings.design_primers(
        seq_args=seq_args,
        global_args=global_args,
    )

    LOGGER.debug(f"Global args: {global_args}")

    LOGGER.debug(f"Primer3 output: {primer3_obj}")

    # Extract primer information
    prim3_res = PrimerSearchResults(primer3_obj=primer3_obj)

    LOGGER.debug(
        "In-silico: context=%s do_insilico_pcr=%s",
        primSet_obj.context,
        getattr(primSet_obj, "do_insilico_pcr", False),
    )
    # run in-silico analysis and populate amplicon summary (optional)
    if isinstance(varInfo_obj, (TranscriptVariantInfo, GenomicVariantInfo)):
        LOGGER.info("Loading primer start and end positions for genomic context")
        prim3_res.load_primer_start_and_end_pos(varInfo_obj)
        if getattr(primSet_obj, "do_insilico_pcr", False):
            LOGGER.info("Running in-silico analysis for designed primer pairs")
            do_insilico_analysis(primSet_obj, prim3_res.primer_pairs)
            LOGGER.info("In-silico analysis completed: ", prim3_res.primer_pairs)
        else:
            LOGGER.info("In-silico PCR disabled; skipping Dicey")
            for pair in prim3_res.primer_pairs:
                pair.amplicons = []
                pair.insilico_status = INSILICO_NOT_APPLICABLE
                pair.insilico_error_detail = None
    elif isinstance(varInfo_obj, SequenceVariantInfo):
        # No genomic coordinates: skip mapped primer positions; Dicey still searches
        # the selected reference when the user enables amplicon check.
        if getattr(primSet_obj, "do_insilico_pcr", False):
            LOGGER.info("Running in-silico analysis for sequence input")
            do_insilico_analysis(primSet_obj, prim3_res.primer_pairs)
        else:
            LOGGER.debug("Sequence input: in-silico PCR disabled; skipping Dicey")
            for pair in prim3_res.primer_pairs:
                pair.amplicons = []
                pair.insilico_status = INSILICO_NOT_APPLICABLE
                pair.insilico_error_detail = None
    else:
        LOGGER.debug(
            "No genomic/transcript variant info; in-silico search not applicable"
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


def _pick_common_reverse_from_primer3(
    raw: dict,
    *,
    min_left_end: int,
    min_gap: int = 80,
) -> tuple[str, int] | None:
    """
    Pick a reverse primer sequence that lies downstream of the discriminating site.

    Primer3 right-primer coordinates are (5' start, length) on the template. For AS-PCR
    the allele-specific forward primer ends at min_left_end; the common reverse must be
    further right so the amplicon spans variant -> reverse.

    min_gap is only a minimum distance to the discriminating site (not the user's full
    product-size minimum, which would reject valid reverses).
    """
    n = int(raw.get("PRIMER_RIGHT_NUM_RETURNED", 0) or 0)
    min_start = int(min_left_end) + int(min_gap)
    for i in range(n):
        pos = raw.get(f"PRIMER_RIGHT_{i}")
        seq = raw.get(f"PRIMER_RIGHT_{i}_SEQUENCE")
        if not pos or not seq:
            continue
        try:
            start = int(pos[0])
        except (TypeError, ValueError, IndexError):
            continue
        if start > min_start:
            return seq, start
    return None


def _run_primer3_with_args(seq_args: dict, global_args: dict) -> PrimerSearchResults:
    if primer3 is None:
        raise RuntimeError("primer3 is not installed in this environment.")
    primer3_obj = primer3.bindings.design_primers(
        seq_args=seq_args,
        global_args=global_args,
    )
    return PrimerSearchResults(primer3_obj=primer3_obj)


def _infer_discriminating_left_ends(var_info: AllelicVariantInfo) -> tuple[int, int]:
    """
    Return (wt_end, mut_end) 0-based coordinates where the left primer must end (3' end).

    SNVs: end at the SNV position.
    Indels: best-effort junction-oriented forcing.
    """
    rel = getattr(var_info, "relative_pos", None)
    if not rel:
        raise ValueError("Allele-specific PCR requires relative_pos on variant info.")
    start = int(rel[0])
    end = int(rel[1])
    indel_type = getattr(var_info, "indel_type", None)
    new_bases = getattr(var_info, "new_bases", "") or ""

    # Default (SNV / unknown): same position in both templates
    wt_end = start
    mut_end = start

    if indel_type == IndelType.INS:
        # For insertions, discrimination depends on the junction and (sometimes) a few
        # downstream bases if the inserted sequence shares a prefix with the WT downstream
        # sequence (e.g. inserted starts with "A" and WT downstream also starts with "A").
        #
        # Find the first position where inserted bases differ from WT downstream bases.
        # Force both primers to end at that discriminating base:
        # - WT ends on the WT base at start+idx (spans junction and reaches mismatch)
        # - MUT ends on the inserted base at start+idx
        wt_downstream = (getattr(var_info, "ref_seq", "") or "")[
            start : start + len(new_bases)
        ]
        mismatch_idx = None
        for i, b in enumerate(new_bases):
            if i >= len(wt_downstream) or b != wt_downstream[i]:
                mismatch_idx = i
                break
        if mismatch_idx is None:
            mismatch_idx = 0

        wt_end = start + mismatch_idx
        mut_end = start + mismatch_idx

    elif indel_type == IndelType.DEL:
        # WT: end at last deleted reference base
        wt_end = end
        # MUT: end immediately before deletion junction
        mut_end = max(0, start - 1)
    elif indel_type == IndelType.DELINS:
        wt_end = end
        mut_end = start + max(0, len(new_bases) - 1)

    return wt_end, mut_end


def primer3_design_allele_specific(
    primSet_obj,
    varInfo_obj: AllelicVariantInfo,
    *,
    max_pairs: int = 5,
) -> dict:
    """
    Design allele-specific PCR primers (WT-specific and MUT-specific).

    Strategy:
    - Use WT template = reference window (varInfo_obj.ref_seq)
    - Use MUT template = varInfo_obj.get_seq('mutated')
    - Pick a common reverse primer on the WT template
    - Run two Primer3 designs with:
        - fixed reverse primer (SEQUENCE_PRIMER_REVCOMP)
        - forced left primer 3' end at discriminating position (SEQUENCE_FORCE_LEFT_END)
        - return primer pairs

    Returns a JSON-serializable dict suitable for persistence.
    """
    if isinstance(varInfo_obj, StructuralVariantWindow):
        raise ValueError(
            "Allele-specific PCR is not supported for structural variants."
        )

    global_args = build_primer3_global_args(primSet_obj)
    global_args = dict(global_args)
    global_args.setdefault("PRIMER_NUM_RETURN", int(max_pairs))

    # WT and MUT templates
    wt_template = getattr(varInfo_obj, "ref_seq", "") or ""
    mut_template = varInfo_obj.get_seq("mutated")
    if not wt_template:
        # Fallback for SequenceVariantInfo etc.
        wt_template = varInfo_obj.get_seq("input").replace("[", "").replace("]", "")

    # Where to force the allele-specific 3' end (left primer end)
    wt_left_end, mut_left_end = _infer_discriminating_left_ends(varInfo_obj)

    # Pick common reverse primer on WT template (no forced left primer)
    common_rev_global = dict(global_args)
    common_rev_global["PRIMER_PICK_LEFT_PRIMER"] = 0
    common_rev_global["PRIMER_PICK_RIGHT_PRIMER"] = 1
    common_rev_global["PRIMER_PICK_INTERNAL_OLIGO"] = 0

    common_rev_seq_args = {
        "SEQUENCE_ID": "as_pcr_common_rev",
        "SEQUENCE_TEMPLATE": wt_template,
    }
    # IMPORTANT: In AS-PCR, the allele-specific primer must overlap the variant.
    # Primer3's SEQUENCE_TARGET excludes primer binding sites in that region, so we
    # must NOT pass SEQUENCE_TARGET here (or for the allele-specific designs).

    # Restrict reverse-primer search to downstream of the discriminating site (WT coords).
    disc_end_wt = wt_left_end
    downstream_start = disc_end_wt + 1
    downstream_len = max(1, len(wt_template) - downstream_start)
    common_rev_seq_args["SEQUENCE_INCLUDED_REGION"] = [downstream_start, downstream_len]

    common_rev_res = _run_primer3_with_args(common_rev_seq_args, common_rev_global)
    raw_common = getattr(common_rev_res, "primer3_obj", None) or {}
    picked = _pick_common_reverse_from_primer3(
        raw_common,
        min_left_end=disc_end_wt,
        min_gap=80,
    )
    common_reverse_seq = None
    common_reverse_start = None
    if picked:
        common_reverse_seq, common_reverse_start = picked
    elif common_rev_res.primer_pairs:
        pair = common_rev_res.primer_pairs[0]
        if (
            pair.right_relPos_start is not None
            and pair.right_relPos_start > disc_end_wt + 80
        ):
            common_reverse_seq = pair.right_seq
            common_reverse_start = pair.right_relPos_end

    if not common_reverse_seq:
        raise ValueError("AS-PCR: could not design a common reverse primer.")

    template_shift = len(mut_template) - len(wt_template)

    def design_for_template(
        template: str, label: str, left_end: int
    ) -> PrimerSearchResults:
        seq_args = {
            "SEQUENCE_ID": f"as_pcr_{label}",
            "SEQUENCE_TEMPLATE": template,
            "SEQUENCE_PRIMER_REVCOMP": common_reverse_seq,
            "SEQUENCE_FORCE_LEFT_END": left_end,
        }
        if (
            label == "mut"
            and template_shift != 0
            and common_reverse_start is not None
            and common_reverse_start > wt_left_end
        ):
            seq_args["SEQUENCE_FORCE_RIGHT_START"] = (
                common_reverse_start + template_shift
            )
        # Do not set SEQUENCE_TARGET in AS-PCR (see note above).

        ga = dict(global_args)
        ga["PRIMER_PICK_LEFT_PRIMER"] = 1
        ga["PRIMER_PICK_RIGHT_PRIMER"] = 1
        ga["PRIMER_PICK_INTERNAL_OLIGO"] = 0
        # With SEQUENCE_FORCE_LEFT_END and a tight size range (18-22), Primer3 only
        # "considers" one candidate per size (5 total). Widen to allow finding a valid
        # primer that still ends exactly at the discriminating base.
        ga["PRIMER_MAX_SIZE"] = max(int(ga.get("PRIMER_MAX_SIZE", 22)), 30)
        # AS-PCR often needs to accept slightly less "perfect" end composition at the
        # forced discriminating position; otherwise Primer3 can reject all candidates.
        ga["PRIMER_GC_CLAMP"] = 0
        ga["PRIMER_MAX_POLY_X"] = max(int(ga.get("PRIMER_MAX_POLY_X", 4)), 6)
        # Widen Tm window a bit for forced-end designs
        if ga.get("PRIMER_OPT_TM") is not None:
            opt_tm = float(ga["PRIMER_OPT_TM"])
            ga["PRIMER_MIN_TM"] = min(
                float(ga.get("PRIMER_MIN_TM", opt_tm - 2)), opt_tm - 5
            )
            ga["PRIMER_MAX_TM"] = max(
                float(ga.get("PRIMER_MAX_TM", opt_tm + 2)), opt_tm + 5
            )
        # Product size can easily fall outside the default 400-800 with a fixed reverse primer.
        # Relax range for AS-PCR to avoid rejecting otherwise valid allele-discriminating primers.
        psr = ga.get("PRIMER_PRODUCT_SIZE_RANGE")
        try:
            if (
                isinstance(psr, list)
                and len(psr) == 2
                and all(isinstance(x, (int, float)) for x in psr)
            ):
                # Extend both ends; transcript/sequence contexts often yield longer products.
                ga["PRIMER_PRODUCT_SIZE_RANGE"] = [
                    min(80, int(psr[0])),
                    max(int(psr[1]), 2000),
                ]
        except Exception:
            pass

        # Sequence input often has higher Tm/hairpin constraints at forced ends; relax mildly for AS-PCR.
        if ga.get("PRIMER_OPT_TM") is not None:
            opt_tm = float(ga["PRIMER_OPT_TM"])
            ga["PRIMER_MAX_TM"] = max(
                float(ga.get("PRIMER_MAX_TM", opt_tm + 5)), opt_tm + 10
            )
        ga.setdefault("PRIMER_MAX_HAIRPIN_TH", 60.0)
        # MUT template at the insertion site often fails strict complementarity checks.
        if label == "mut":
            ga["PRIMER_MAX_SELF_ANY"] = max(int(ga.get("PRIMER_MAX_SELF_ANY", 12)), 20)
            ga["PRIMER_MAX_SELF_END"] = max(int(ga.get("PRIMER_MAX_SELF_END", 8)), 16)
            ga["PRIMER_PAIR_MAX_COMPL_ANY"] = max(
                int(ga.get("PRIMER_PAIR_MAX_COMPL_ANY", 12)), 20
            )
            ga["PRIMER_PAIR_MAX_COMPL_END"] = max(
                int(ga.get("PRIMER_PAIR_MAX_COMPL_END", 8)), 16
            )
            ga["PRIMER_PICK_ANYWAY"] = 1
        return _run_primer3_with_args(seq_args, ga)

    wt_res = design_for_template(wt_template, "wt", wt_left_end)
    mut_res = design_for_template(mut_template, "mut", mut_left_end)

    # Populate in-silico + mapped positions like standard design
    for res in (wt_res, mut_res):
        if isinstance(varInfo_obj, (TranscriptVariantInfo, GenomicVariantInfo)):
            res.load_primer_start_and_end_pos(varInfo_obj)
            if getattr(primSet_obj, "do_insilico_pcr", False):
                do_insilico_analysis(primSet_obj, res.primer_pairs)
            else:
                for pair in res.primer_pairs:
                    pair.amplicons = []
                    pair.insilico_status = INSILICO_NOT_APPLICABLE
                    pair.insilico_error_detail = None
        elif isinstance(varInfo_obj, SequenceVariantInfo):
            if getattr(primSet_obj, "do_insilico_pcr", False):
                do_insilico_analysis(primSet_obj, res.primer_pairs)
            else:
                for pair in res.primer_pairs:
                    pair.amplicons = []
                    pair.insilico_status = INSILICO_NOT_APPLICABLE
                    pair.insilico_error_detail = None
        else:
            for pair in res.primer_pairs:
                pair.amplicons = []
                pair.insilico_status = INSILICO_NOT_APPLICABLE
                pair.insilico_error_detail = None

    return {
        "design_type": "allele_specific",
        "common_reverse_primer": common_reverse_seq,
        "wt_left_force_end": wt_left_end,
        "mut_left_force_end": mut_left_end,
        "wt": wt_res.to_dict(),
        "mut": mut_res.to_dict(),
    }
