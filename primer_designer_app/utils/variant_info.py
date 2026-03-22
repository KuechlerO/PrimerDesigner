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
        """Parse variant information from sequence annotation"""
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
