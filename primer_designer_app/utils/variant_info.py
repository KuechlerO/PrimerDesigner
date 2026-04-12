from dataclasses import dataclass
from typing import Optional, Tuple
import logging

from enum import Enum
from primer_designer_app.utils.ensembl_client import EnsemblClient


LOGGER = logging.getLogger(__name__)


VARIANT_FLANKING = 1000


class IndelType(Enum):
    INS = 'INS'
    DEL = 'DEL'
    DELINS = 'DelIns'
    SNV = 'SNV'
    NONE = ''


class ReferenceType(Enum):
    CDS = 'cds'
    CDNA = 'cdna'
    PROTEIN = 'protein'
    NONE = ''


@dataclass
class VariantInfo:
    """Information for a variant (provided by user throught web app interface)"""

    ref_seq: str = ''
    ref_bases: str = ''
    new_bases: str = ''
    gene_ID: str = ''
    gene_symbol: str = ''
    ref_genome: str = ''
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
            field
        ) in self.__dataclass_fields__:  # Automatically get all dataclass fields
            if field in kwargs:
                setattr(self, field, kwargs[field])
        # Ignore unknown fields (e.g., AS_seq from TranscriptVariantInfo)
        LOGGER.debug(
            f"Ignored unexpected kwargs: {set(kwargs) - set(self.__dataclass_fields__.keys())}"
        )
        self.__post_init__()

    def __post_init__(self):
        if self.relative_pos is None:
            raise ValueError('relative_pos must be provided for VariantInfo')

    def set_attribute(self, key, value):
        setattr(self, key, value)

    def _determine_indel_type(self):
        ref_length = len(self.ref_bases) if self.ref_bases.isalpha() else 0
        new_length = len(self.new_bases) if self.new_bases.isalpha() else 0

        LOGGER.debug(
            f"Determining indel type with ref_bases: '{self.ref_bases}', new_bases: '{self.new_bases}', ref_length: {ref_length}, new_length: {new_length}"
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
                'Could not determine indel type from provided bases and positions.'
            )

    def _load_geneDetails(self, ensembl_client: EnsemblClient) -> Tuple[str, str]:
        """Use ensembl_client to load gene symbols and IDs based on genomic position

        Returns:
            Optional[str]: gene ID and symbol tuple if found, else None
        """

        if not self.genomic_pos or not self.genomic_pos.get('pos'):
            raise ValueError(
                "Genomic position with 'pos' key must be provided to load gene details"
            )

        chr_label = self.genomic_pos['chr']
        if chr_label == '23':
            chr_label = 'X'
        elif chr_label == '24':
            chr_label = 'Y'

        start = self.genomic_pos['pos'][0]
        if len(self.genomic_pos['pos']) > 1:
            end = self.genomic_pos['pos'][-1]
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
            LOGGER.exception('Ensembl lookup failed')
        return ('', '')

    def get_genomic_pos(self):
        return self.genomic_pos['pos']

    def get_seq(self, output_type: str) -> str:
        """Construct the mutated or input sequence based on the reference sequence and variant information

        Args:
            output_type (str): "mutated" for the sequence with the variant applied, "input" for the original input format with brackets

        Returns:
            str: The constructed sequence based on the specified output type
        """
        if output_type == 'mutated':
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
        elif output_type == 'input':
            if self.indel_type == IndelType.SNV:
                return (
                    self.ref_seq[: self.relative_pos[0]]
                    + '['
                    + self.ref_bases
                    + '>'
                    + self.new_bases
                    + ']'
                    + self.ref_seq[self.relative_pos[0] + 1 :]
                )
            elif self.indel_type == IndelType.INS:
                return (
                    self.ref_seq[: self.relative_pos[0]]
                    + '[-/'
                    + self.new_bases
                    + ']'
                    + self.ref_seq[self.relative_pos[0] :]
                )
            elif self.indel_type == IndelType.DEL:
                return (
                    self.ref_seq[: self.relative_pos[0]]
                    + '['
                    + self.ref_bases
                    + '/-]'
                    + self.ref_seq[self.relative_pos[1] + 1 :]
                )
            elif self.indel_type == IndelType.DELINS:
                return (
                    self.ref_seq[: self.relative_pos[0]]
                    + '['
                    + self.ref_bases
                    + '/'
                    + self.new_bases
                    + ']'
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
            f"Running GenomicVariantInfo init with relative_pos: {relative_pos}, genomic_pos: {genomic_pos}"
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
                'Genomic position must be provided to fetch sequence snippet'
            )

        start_position = max(1, self.genomic_pos['pos'][0] - flank)
        end_position = self.genomic_pos['pos'][-1] + flank
        seq = ensembl_client.get_genomic_sequence(
            self.genomic_pos['chr'], start_position, end_position
        )
        return seq


class TranscriptVariantInfo(VariantInfo):
    """Information for a variant when gene or transcript ID is provided"""

    transcript_id: str

    def __init__(self, transcript_id: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.transcript_id = transcript_id

        LOGGER.debug(
            f"2. Relative position: {self.relative_pos}, Reference type: {self.reference_type}"
        )

        # Fetch data from Ensembl
        ensembl_client = EnsemblClient(ref_genome=self.ref_genome)
        self.gene_symbol, used_transcript_id = ensembl_client.get_gene_symbol_for_transcriptID(
            self.transcript_id
        )
        self.transcript_id = used_transcript_id

        self.ref_seq = ensembl_client.get_transcript_sequence(
            self.transcript_id, self.reference_type.value
        )

        LOGGER.debug(
            f"Before extracting ref_bases, ref_seq: {self.ref_seq}, relative_pos: {self.relative_pos}, ref_bases: {self.ref_bases}"
        )
        if not self.ref_bases or self.ref_bases == '':
            self.ref_bases = self.ref_seq[
                self.relative_pos[0] : self.relative_pos[1] + 1
            ]
            LOGGER.debug(
                f"Extracting ref_bases from ref_seq using relative_pos: {self.relative_pos}"
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
            f"Mapping transcript coordinates to genomic coordinates for transcript ID: {self.transcript_id}"
        )

        data = ensembl_client.map_coordinates(
            self.transcript_id,
            self.relative_pos[0],
            self.relative_pos[1],
            self.reference_type.value,
        )

        if len(data['mappings']) == 1:
            pos = [
                int(data['mappings'][0]['start']) - 1,
                int(data['mappings'][0]['end']) - 1,
            ]
        else:
            pos = []
            for i in range(len(data['mappings'])):
                pos.append(
                    [int(data['mappings'][i]['start']), int(data['mappings'][i]['end'])]
                )

        chromosome = data['mappings'][0]['seq_region_name']
        strand_type = 'antisense' if data['mappings'][0]['strand'] == -1 else 'sense'

        genomic_pos = {
            'chr': chromosome,
            'pos': pos,
            'strand_type': strand_type,
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
        input_seq = ''.join(input_seq.split())
        if '[' not in input_seq or ']' not in input_seq:
            raise ValueError('Sequence annotation missing brackets')

        start = input_seq.index('[')
        end = input_seq.index(']')
        var_input = input_seq[start + 1 : end]

        if '>' in var_input:  # SNV like [A>G]
            left, right = var_input.split('>')
            ref_bases = left
            new_bases = right
            relative_pos = (start, start)
        elif '/' in var_input:  # indel like [A/T] or [-/T]
            left, right = var_input.split('/')
            ref_bases = left if left not in ['-', None, ''] else ''
            new_bases = right if right not in ['-', None, ''] else ''

            LOGGER.debug(
                f"Calculating relative position for deletion with ref_bases: {ref_bases}, new_bases: {new_bases}"
            )
            if ref_bases == '' and new_bases != '':
                relative_pos = (start, start)  # Insertion at this position
            else:
                # Deletion of these bases
                relative_pos = (start, start + len(ref_bases) - 1)
        else:
            raise ValueError(f"Invalid variant annotation format: {var_input}")

        ref_seq = input_seq.replace(f"[{var_input}]", ref_bases)
        return ref_seq, ref_bases, new_bases, relative_pos





@dataclass
class StructuralVariantInfo:
    """
    Repräsentiert entweder
    1) die gesamte Strukturvariante oder
    2) ein konkretes Designfenster dieser Strukturvariante.

    Die Klasse ist bewusst unabhängig von VariantInfo.
    """

    chromosome: str
    start_position: int
    end_position: int
    structural_variant_type: str
    reference_genome: str = 'GRCh37'

    # Fensterspezifische Informationen
    label: str = 'structural_variant'
    window_start_genomic: Optional[int] = None
    window_end_genomic: Optional[int] = None

    # Primer3-relevante Informationen
    window_sequence: str = ''
    target_start_in_window: Optional[int] = None
    target_end_in_window: Optional[int] = None
    target_length: Optional[int] = None

    def __post_init__(self) -> None:
        self.chromosome = self._normalize_chromosome(self.chromosome)
        self.structural_variant_type = self.structural_variant_type.strip().lower()

        if self.start_position < 1 or self.end_position < 1:
            raise ValueError('Start and end position must be positive integers')

        if self.start_position > self.end_position:
            raise ValueError('Start position must not be greater than end position')

        allowed_sv_types = {'deletion', 'duplication'}
        if self.structural_variant_type not in allowed_sv_types:
            raise ValueError(
                f"Unsupported structural variant type: {self.structural_variant_type}"
            )

        if self.window_start_genomic is None:
            self.window_start_genomic = self.start_position

        if self.window_end_genomic is None:
            self.window_end_genomic = self.end_position

        if self.window_start_genomic < 1:
            self.window_start_genomic = 1

        if self.window_start_genomic > self.window_end_genomic:
            raise ValueError('Window start must not be greater than window end')

    @staticmethod
    def _normalize_chromosome(chromosome: str) -> str:
        chromosome = chromosome.strip().upper()
        if chromosome.startswith('CHR'):
            chromosome = chromosome[3:]
        return chromosome

    @property
    def structural_variant_length(self) -> int:
        return self.end_position - self.start_position + 1

    def get_genomic_pos(self) -> list[int]:
        """
        Kompatible Methode für bestehenden Code.
        Gibt die genomischen Koordinaten des aktuellen Designfensters zurück.
        """
        return [self.window_start_genomic, self.window_end_genomic]

    def get_seq(self, output_type: str) -> str:
        """
        Kompatible Methode für primer3_design_primers(...).

        Für SV gibt es in dieser ersten Version nur die Fenstersequenz,
        daher unterstützen wir nur 'mutated'.
        """
        if output_type != 'mutated':
            raise ValueError(
                "StructuralVariantInfo supports only output_type='mutated'"
            )

        if not self.window_sequence:
            raise ValueError('window_sequence has not been loaded yet')

        return self.window_sequence

    def load_window_sequence(self, mask_feature: str = '0', strand: str = '1') -> str:
        """
        Holt die Referenzsequenz für das aktuelle Fenster von Ensembl
        und speichert sie in window_sequence.
        """
        ensembl_client = EnsemblClient(ref_genome=self.reference_genome)
        self.window_sequence = ensembl_client.get_genomic_sequence(
            chromosome=self.chromosome,
            start=self.window_start_genomic,
            end=self.window_end_genomic,
            strand=strand,
            mask_feature=mask_feature,
        )

        return self.window_sequence

    def define_target_in_window(
        self,
        requested_target_length: int = 200,
        minimum_margin: int = 10,
    ) -> list[int]:
        """
        Definiert einen Zielbereich innerhalb des geladenen Fensters.

        Wichtig:
        - intern speichern wir Start und Ende
        - weil primer_settings.set_target(...) im Projekt offenbar
        mit (start, end) arbeitet
        """
        if not self.window_sequence:
            raise ValueError('window_sequence must be loaded before defining a target')

        sequence_length = len(self.window_sequence)
        max_possible_target_length = sequence_length - (2 * minimum_margin)

        if max_possible_target_length <= 0:
            raise ValueError(
                f"Window sequence too short for primer design: length={sequence_length}"
            )

        self.target_length = min(requested_target_length, max_possible_target_length)

        if self.target_length <= 0:
            raise ValueError('Target length must be positive')

        self.target_start_in_window = (sequence_length - self.target_length) // 2
        self.target_end_in_window = (
            self.target_start_in_window + self.target_length - 1
        )

        return [self.target_start_in_window, self.target_end_in_window]

    def get_target_interval_in_window(self) -> list[int]:
        if self.target_start_in_window is None or self.target_end_in_window is None:
            raise ValueError('Target positions have not been defined yet')

        if self.target_end_in_window < self.target_start_in_window:
            raise ValueError('Target end must not be smaller than target start')

        return [self.target_start_in_window, self.target_end_in_window]


    def _center_window_in_interval(
        self,
        interval_start: int,
        interval_end: int,
        window_size: int,
    ) -> tuple[int, int]:
        """
        Platziert ein Fenster der Länge window_size mittig in einem Intervall.
        """
        interval_length = interval_end - interval_start + 1

        if window_size > interval_length:
            raise ValueError(
                f"Window size {window_size} is larger than interval length {interval_length}"
            )

        window_start = interval_start + (interval_length - window_size) // 2
        window_end = window_start + window_size - 1
        return window_start, window_end


    def create_design_windows(
        self,
        flank_window_size: int = 5000,
        internal_window_size: int = 2000,
        minimum_internal_window_size: int = 120,
    ) -> list['StructuralVariantInfo']:
        """
        Erzeugt vier Designfenster:
        - upstream
        - downstream
        - internal_1
        - internal_2

        Die beiden internen Fenster werden bewusst in die linke und rechte
        Hälfte der SV gelegt, damit sie sich nicht überschneiden.
        """

        sv_length = self.structural_variant_length

        # --- Äußere Fenster ---
        upstream_start = max(1, self.start_position - flank_window_size)
        upstream_end = self.start_position - 1

        downstream_start = self.end_position + 1
        downstream_end = self.end_position + flank_window_size

        if upstream_start > upstream_end:
            raise ValueError('No valid upstream window can be created')

        if downstream_start > downstream_end:
            raise ValueError('No valid downstream window can be created')

        upstream_window = StructuralVariantInfo(
            chromosome=self.chromosome,
            start_position=self.start_position,
            end_position=self.end_position,
            structural_variant_type=self.structural_variant_type,
            reference_genome=self.reference_genome,
            label='upstream',
            window_start_genomic=upstream_start,
            window_end_genomic=upstream_end,
        )

        downstream_window = StructuralVariantInfo(
            chromosome=self.chromosome,
            start_position=self.start_position,
            end_position=self.end_position,
            structural_variant_type=self.structural_variant_type,
            reference_genome=self.reference_genome,
            label='downstream',
            window_start_genomic=downstream_start,
            window_end_genomic=downstream_end,
        )

        # --- Innere Fenster: linke und rechte Hälfte der SV ---
        left_half_start = self.start_position
        left_half_end = self.start_position + (sv_length // 2) - 1

        right_half_start = left_half_end + 1
        right_half_end = self.end_position

        left_half_length = left_half_end - left_half_start + 1
        right_half_length = right_half_end - right_half_start + 1

        # Jedes interne Fenster muss in seine Hälfte passen.
        effective_internal_window_size = min(
            internal_window_size,
            left_half_length,
            right_half_length,
        )

        # Wenn die SV zu klein ist, würden die internen Fenster zwar formal passen,
        # aber biologisch/technisch oft keinen Sinn mehr ergeben.
        if effective_internal_window_size < minimum_internal_window_size:
            raise ValueError(
                'Structural variant is too short for two meaningful non-overlapping '
                f"internal windows. Maximum possible size per internal window is "
                f"{effective_internal_window_size} bp."
            )

        internal_1_start, internal_1_end = self._center_window_in_interval(
            interval_start=left_half_start,
            interval_end=left_half_end,
            window_size=effective_internal_window_size,
        )

        internal_2_start, internal_2_end = self._center_window_in_interval(
            interval_start=right_half_start,
            interval_end=right_half_end,
            window_size=effective_internal_window_size,
        )

        internal_window_1 = StructuralVariantInfo(
            chromosome=self.chromosome,
            start_position=self.start_position,
            end_position=self.end_position,
            structural_variant_type=self.structural_variant_type,
            reference_genome=self.reference_genome,
            label='internal_1',
            window_start_genomic=internal_1_start,
            window_end_genomic=internal_1_end,
        )

        internal_window_2 = StructuralVariantInfo(
            chromosome=self.chromosome,
            start_position=self.start_position,
            end_position=self.end_position,
            structural_variant_type=self.structural_variant_type,
            reference_genome=self.reference_genome,
            label='internal_2',
            window_start_genomic=internal_2_start,
            window_end_genomic=internal_2_end,
        )

        return [
            upstream_window,
            downstream_window,
            internal_window_1,
            internal_window_2,
        ]

    def prepare_for_primer_design(
        self,
        requested_target_length: int = 200,
        minimum_margin: int = 10,
        mask_feature: str = '0',
        strand: str = '1',
    ) -> None:
        """
        Komfortmethode:
        1. Sequenz laden
        2. Primer3-Zielbereich definieren
        """
        self.load_window_sequence(mask_feature=mask_feature, strand=strand)
        self.define_target_in_window(
            requested_target_length=requested_target_length,
            minimum_margin=minimum_margin,
        )
