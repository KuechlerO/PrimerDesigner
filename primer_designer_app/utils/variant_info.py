from dataclasses import dataclass, field
from typing import Optional, Tuple, List
import logging

from enum import Enum
from primer_designer_app.utils.ensembl_client import EnsemblClient

LOGGER = logging.getLogger(__name__)


VARIANT_FLANKING = 1000


class IndelType(Enum):
    INS = "INS"
    DEL = "DEL"
    DELINS = "DelIns"
    SNV = "SNV"
    NONE = ""


class ReferenceType(Enum):
    CDS = "cds"
    CDNA = "cdna"
    PROTEIN = "protein"
    NONE = ""


@dataclass
class VariantInfo:
    """Information for a variant (provided by user throught web app interface)"""

    ref_seq: str = ""
    ref_bases: str = ""
    new_bases: str = ""
    gene_ID: str = ""
    gene_symbol: str = ""
    ref_genome: str = ""
    indel_type: IndelType = IndelType.NONE
    # relative position within seq (0-based, inclusive)
    relative_pos: Optional[Tuple[int, int]] = (
        None  # (start, end) positions of the variant within the sequence
    )
    genomic_pos: Optional[dict] = None
    reference_type: ReferenceType = ReferenceType.NONE

    def __init__(self, **kwargs):
        # Set attributes for known fields
        for (
            field_name
        ) in self.__dataclass_fields__:  # Automatically get all dataclass fields
            if field_name in kwargs:
                setattr(self, field_name, kwargs[field_name])
        # Ignore unknown fields (e.g., AS_seq from TranscriptVariantInfo)
        LOGGER.debug(
            f"Ignored unexpected kwargs: "
            f"{set(kwargs) - set(self.__dataclass_fields__.keys())}"
        )
        self.__post_init__()

    def __post_init__(self):
        if self.relative_pos is None:
            raise ValueError("relative_pos must be provided for VariantInfo")

    def set_attribute(self, key, value):
        setattr(self, key, value)

    def _determine_indel_type(self):
        ref_length = len(self.ref_bases) if self.ref_bases.isalpha() else 0
        new_length = len(self.new_bases) if self.new_bases.isalpha() else 0

        LOGGER.debug(
            f"Determining indel type with ref_bases: '{self.ref_bases}', "
            f"new_bases: '{self.new_bases}', ref_length: {ref_length}, "
            f"new_length: {new_length}"
        )

        if ref_length == 0 and new_length >= 1:
            return IndelType.INS
        elif ref_length == 1 and new_length == 1:
            return IndelType.SNV
        elif ref_length >= 1:
            if new_length >= 1:
                return IndelType.DELINS
            elif not new_length or new_length == 0:
                return IndelType.DEL
        else:
            raise ValueError(
                "Could not determine indel type from provided bases and positions."
            )

    def _load_geneDetails(self, ensembl_client: EnsemblClient) -> Tuple[str, str]:
        """Use ensembl_client to load gene symbols and IDs based on genomic position

        Returns:
            Optional[str]: gene ID and symbol tuple if found, else None
        """

        if not self.genomic_pos or not self.genomic_pos.get("pos"):
            raise ValueError(
                "Genomic position with 'pos' key must be provided to load gene details"
            )

        chr_label = self.genomic_pos["chr"]
        if chr_label == "23":
            chr_label = "X"
        elif chr_label == "24":
            chr_label = "Y"

        start = self.genomic_pos["pos"][0]
        if len(self.genomic_pos["pos"]) > 1:
            end = self.genomic_pos["pos"][-1]
        else:
            end = start

        try:
            genes = ensembl_client.get_overlapped_genes_details_for_region(
                chr_label, start, end
            )
            if genes:
                # Only return the first gene found
                for gene_ID, gene_symbol in genes.items():
                    if gene_ID and gene_symbol:
                        return (gene_ID, gene_symbol)
        except Exception:
            LOGGER.exception("Ensembl lookup failed")
        return ("", "")

    def get_genomic_pos(self):
        return self.genomic_pos["pos"]

    def get_seq(self, output_type: str) -> str:
        """Construct the mutated or input sequence based on the reference
        sequence and variant information

        Args:
            output_type (str): "mutated" for the sequence with the variant
            applied, "input" for the original input format with brackets

        Returns:
            str: The constructed sequence based on the specified output type
        """
        if output_type == "mutated":
            if self.indel_type == IndelType.SNV:
                return (
                    self.ref_seq[: self.relative_pos[0]]
                    + self.new_bases
                    + self.ref_seq[self.relative_pos[0] + 1 :]
                )
            elif self.indel_type == IndelType.INS:
                return (
                    self.ref_seq[: self.relative_pos[0]]
                    + self.new_bases
                    + self.ref_seq[self.relative_pos[0] :]
                )
            elif self.indel_type == IndelType.DEL:
                return (
                    self.ref_seq[: self.relative_pos[0]]
                    + self.ref_seq[self.relative_pos[1] + 1 :]
                )
            elif self.indel_type == IndelType.DELINS:
                return (
                    self.ref_seq[: self.relative_pos[0]]
                    + self.new_bases
                    + self.ref_seq[self.relative_pos[1] + 1 :]
                )
            else:
                raise ValueError(f"Invalid indel type: {self.indel_type}")
        elif output_type == "input":
            if self.indel_type == IndelType.SNV:
                return (
                    self.ref_seq[: self.relative_pos[0]]
                    + "["
                    + self.ref_bases
                    + ">"
                    + self.new_bases
                    + "]"
                    + self.ref_seq[self.relative_pos[0] + 1 :]
                )
            elif self.indel_type == IndelType.INS:
                return (
                    self.ref_seq[: self.relative_pos[0]]
                    + "[-/"
                    + self.new_bases
                    + "]"
                    + self.ref_seq[self.relative_pos[0] :]
                )
            elif self.indel_type == IndelType.DEL:
                return (
                    self.ref_seq[: self.relative_pos[0]]
                    + "["
                    + self.ref_bases
                    + "/-]"
                    + self.ref_seq[self.relative_pos[1] + 1 :]
                )
            elif self.indel_type == IndelType.DELINS:
                return (
                    self.ref_seq[: self.relative_pos[0]]
                    + "["
                    + self.ref_bases
                    + "/"
                    + self.new_bases
                    + "]"
                    + self.ref_seq[self.relative_pos[1] + 1 :]
                )
            else:
                raise ValueError(f"Invalid indel type: {self.indel_type}")
        else:
            raise ValueError(f"Invalid output type: {output_type}")


class GenomicVariantInfo(VariantInfo):
    def __init__(
        self,
        relative_pos: Optional[Tuple[int, int]] = [VARIANT_FLANKING, VARIANT_FLANKING],
        genomic_pos: Optional[dict] = None,
        *args,
        **kwargs,
    ):
        LOGGER.info(
            f"Running GenomicVariantInfo init with relative_pos: "
            f"{relative_pos}, genomic_pos: {genomic_pos}"
        )
        # Pass only the arguments expected by the parent class
        super().__init__(
            relative_pos=relative_pos, genomic_pos=genomic_pos, *args, **kwargs
        )
        ensembl_client = EnsemblClient(ref_genome=self.ref_genome)

        # Load sequence snippet if genomic position is provided
        self.ref_seq = self._get_sequence_snippet(ensembl_client)
        self.ref_bases = self.ref_seq[self.relative_pos[0] : self.relative_pos[1] + 1]

        # Load gene details if genomic position is provided
        if self.genomic_pos and not self.gene_ID:
            self.gene_ID, self.gene_symbol = self._load_geneDetails(ensembl_client)

        # Determine indel type if not set
        self.indel_type = self._determine_indel_type()

    def _get_sequence_snippet(
        self, ensembl_client: EnsemblClient, flank: int = VARIANT_FLANKING
    ) -> str:
        """Fetch seq snippet using a Django model (injected) with a substring helper."""
        if not self.genomic_pos:
            raise ValueError(
                "Genomic position must be provided to fetch sequence snippet"
            )

        start_position = max(1, self.genomic_pos["pos"][0] - flank)
        end_position = self.genomic_pos["pos"][-1] + flank
        seq = ensembl_client.get_genomic_sequence(
            self.genomic_pos["chr"], start_position, end_position
        )
        return seq


class TranscriptVariantInfo(VariantInfo):
    """Information for a variant when gene or transcript ID is provided"""

    transcript_id: str

    def __init__(self, transcript_id: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.transcript_id = transcript_id

        LOGGER.debug(
            f"2. Relative position: {self.relative_pos}, "
            f"Reference type: {self.reference_type}"
        )

        # Fetch data from Ensembl
        ensembl_client = EnsemblClient(ref_genome=self.ref_genome)
        self.gene_symbol, used_transcript_id = (
            ensembl_client.get_gene_symbol_for_transcriptID(self.transcript_id)
        )
        self.transcript_id = used_transcript_id

        self.ref_seq = ensembl_client.get_transcript_sequence(
            self.transcript_id, self.reference_type.value
        )

        LOGGER.debug(
            f"Before extracting ref_bases, ref_seq: {self.ref_seq}, "
            f"relative_pos: {self.relative_pos}, ref_bases: {self.ref_bases}"
        )
        if not self.ref_bases or self.ref_bases == "":
            self.ref_bases = self.ref_seq[
                self.relative_pos[0] : self.relative_pos[1] + 1
            ]
            LOGGER.debug(
                f"Extracting ref_bases from ref_seq using relative_pos: "
                f"{self.relative_pos}"
            )
            LOGGER.debug(
                f"Extracted ref_bases: {self.ref_bases} from ref_seq: {self.ref_seq}"
            )
        self.genomic_pos = self._get_genomic_pos(ensembl_client)

        # Load gene details if genomic position is provided
        if self.genomic_pos and not self.gene_ID:
            self.gene_ID, self.gene_symbol = self._load_geneDetails(ensembl_client)

        # Determine indel type if not set
        self.indel_type = self._determine_indel_type()

    def _get_genomic_pos(self, ensembl_client: EnsemblClient) -> dict:
        """Get genomic position given the transcript ID and relative positions

        Returns:
            genomic_pos (dict): chromosome, position, strand_type
        """
        LOGGER.info(f"self.relative_pos: {self.relative_pos}")
        LOGGER.info(f"self.reference_type: {self.reference_type}")
        LOGGER.info(
            f"Mapping transcript coordinates to genomic coordinates "
            f"for transcript ID: {self.transcript_id}"
        )

        data = ensembl_client.map_coordinates(
            self.transcript_id,
            self.relative_pos[0],
            self.relative_pos[1],
            self.reference_type.value,
        )

        if len(data["mappings"]) == 1:
            pos = [
                int(data["mappings"][0]["start"]) - 1,
                int(data["mappings"][0]["end"]) - 1,
            ]
        else:
            pos = []
            for i in range(len(data["mappings"])):
                pos.append(
                    [int(data["mappings"][i]["start"]), int(data["mappings"][i]["end"])]
                )

        chromosome = data["mappings"][0]["seq_region_name"]
        strand_type = "antisense" if data["mappings"][0]["strand"] == -1 else "sense"

        genomic_pos = {
            "chr": chromosome,
            "pos": pos,
            "strand_type": strand_type,
        }
        return genomic_pos


class SequenceVariantInfo(VariantInfo):
    """Parse inline sequence annotations, e.g. AC[3>A]GT or ACG[-/T]T"""

    def __init__(self, input_seq, **kwargs):
        self.ref_seq, self.ref_bases, self.new_bases, self.relative_pos = (
            self._parse_input_sequence(input_seq)
        )
        LOGGER.debug(f"3. Relative position: {self.relative_pos}")
        LOGGER.debug(f"Variant info before: {self}")
        super().__init__(**kwargs)
        self.indel_type = self._determine_indel_type()
        LOGGER.debug(f"Variant info after: {self}")

    def _parse_input_sequence(self, input_seq: str):
        """Parse variant information from sequence annotation.

        Strips all whitespace (spaces, line breaks, tabs) so pasted multi-line
        sequences are accepted.
        """
        input_seq = "".join(input_seq.split())
        if "[" not in input_seq or "]" not in input_seq:
            raise ValueError("Sequence annotation missing brackets")

        start = input_seq.index("[")
        end = input_seq.index("]")
        var_input = input_seq[start + 1 : end]

        if ">" in var_input:  # SNV like [A>G]
            left, right = var_input.split(">")
            ref_bases = left
            new_bases = right
            relative_pos = (start, start)
        elif "/" in var_input:  # indel like [A/T] or [-/T]
            left, right = var_input.split("/")
            ref_bases = left if left not in ["-", None, ""] else ""
            new_bases = right if right not in ["-", None, ""] else ""

            LOGGER.debug(
                f"Calculating relative position for deletion with ref_bases: "
                f"{ref_bases}, new_bases: {new_bases}"
            )
            if ref_bases == "" and new_bases != "":
                relative_pos = (start, start)  # Insertion at this position
            else:
                # Deletion of these bases
                relative_pos = (start, start + len(ref_bases) - 1)
        else:
            raise ValueError(f"Invalid variant annotation format: {var_input}")

        ref_seq = input_seq.replace(f"[{var_input}]", ref_bases)
        return ref_seq, ref_bases, new_bases, relative_pos


@dataclass
class StructuralVariantWindow:
    label: str
    window_start_genomic: int
    window_end_genomic: int
    target_start_in_window: Optional[int] = None
    target_length: Optional[int] = None
    window_sequence: str = ""

    @property
    def window_length(self) -> int:
        return self.window_end_genomic - self.window_start_genomic + 1

    def set_target(self, target_start_in_window: int, target_length: int) -> None:
        if target_start_in_window < 0:
            raise ValueError("target_start_in_window must be >= 0")
        if target_length <= 0:
            raise ValueError("target_length must be > 0")
        if target_start_in_window + target_length > self.window_length:
            raise ValueError("Target exceeds window boundaries")

        self.target_start_in_window = target_start_in_window
        self.target_length = target_length

    def set_default_target(self, target_length: int = 150) -> None:
        effective_target_length = min(target_length, self.window_length)
        target_start_in_window = (self.window_length - effective_target_length) // 2
        self.set_target(target_start_in_window, effective_target_length)

    def get_primer3_target(self) -> list[int]:
        if self.target_start_in_window is None or self.target_length is None:
            raise ValueError("Primer3 target has not been set")
        return [self.target_start_in_window, self.target_length]

    def load_window_sequence(self, chromosome: str, reference_genome: str) -> None:
        client = EnsemblClient(ref_genome=reference_genome)
        self.window_sequence = client.get_genomic_sequence(
            chromosome=chromosome,
            start=self.window_start_genomic,
            end=self.window_end_genomic,
            mask_feature="0",
        )

    def get_seq(self, output_type: str) -> str:
        """
        Return the loaded window sequence for Primer3 compatibility.
        """
        if output_type != "mutated":
            raise ValueError(
                "StructuralVariantWindow supports only output_type='mutated'"
            )

        if not self.window_sequence:
            raise ValueError("window_sequence has not been loaded yet")

        return self.window_sequence

    def get_genomic_pos(self) -> list[int]:
        return [self.window_start_genomic, self.window_end_genomic]


@dataclass
class StructuralVariantInfo:
    chromosome: str
    start_position: int
    end_position: int
    reference_genome: str
    windows: List[StructuralVariantWindow] = field(default_factory=list)

    def __post_init__(self) -> None:
        # Validate core coordinate invariants.
        if self.start_position <= 0 or self.end_position <= 0:
            raise ValueError("Genomic positions must be positive integers")

        if self.start_position > self.end_position:
            raise ValueError("Start position must not exceed end position")

    @property
    def structural_variant_length(self) -> int:
        return self.end_position - self.start_position + 1

    def create_design_windows(
        self,
        flank_window_size: int = 5000,
        internal_window_size: int = 2000,
    ) -> List[StructuralVariantWindow]:
        sv_length = self.structural_variant_length

        if sv_length < 50:
            raise ValueError("Structural variant must span at least 50 bases")

        upstream_window_end = self.start_position - 1
        if upstream_window_end < 1:
            raise ValueError("Upstream window cannot be created at chromosome start")

        half_sv_length = sv_length // 2
        internal_size = min(internal_window_size, half_sv_length)

        if internal_size < 25:
            raise ValueError("Internal window size is too small for primer design")

        upstream_window = StructuralVariantWindow(
            label="upstream",
            window_start_genomic=max(1, self.start_position - flank_window_size),
            window_end_genomic=upstream_window_end,
        )

        downstream_window = StructuralVariantWindow(
            label="downstream",
            window_start_genomic=self.end_position + 1,
            window_end_genomic=self.end_position + flank_window_size,
        )

        internal_1 = StructuralVariantWindow(
            label="internal_1",
            window_start_genomic=self.start_position,
            window_end_genomic=self.start_position + internal_size - 1,
        )

        internal_2 = StructuralVariantWindow(
            label="internal_2",
            window_start_genomic=self.end_position - internal_size + 1,
            window_end_genomic=self.end_position,
        )

        self.windows = [
            upstream_window,
            downstream_window,
            internal_1,
            internal_2,
        ]
        return self.windows
