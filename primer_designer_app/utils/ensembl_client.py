import requests
from requests.adapters import HTTPAdapter, Retry
from primer_designer_app.exceptions import (
    InvalidTranscriptIdError,
    InvalidTranscriptVersionError,
)

server37 = 'https://grch37.rest.ensembl.org'
server38 = 'http://rest.ensembl.org'

import logging

logger = logging.getLogger(__name__)


class EnsemblClient:
    def __init__(self, ref_genome: str = 'GRCh38'):
        if ref_genome == 'GRCh38':
            self.server = server38
        elif ref_genome == 'GRCh37':
            self.server = server37
        else:
            raise ValueError(
                f"Unsupported reference genome: {ref_genome}. Supported values are 'GRCh38' and 'GRCh37'."
            )

        self.session = requests.Session()
        retries = Retry(
            total=3, backoff_factor=0.3, status_forcelist=[429, 500, 502, 503, 504]
        )
        self.session.mount('https://', HTTPAdapter(max_retries=retries))

    def split_transcript_id(self, transcript_id: str):
        transcript_id = transcript_id.strip()

        if '.' in transcript_id:
            base_id, version = transcript_id.split('.', 1)
            return base_id, int(version)

        return transcript_id, None

    def get_transcript_sequence(
        self, transcript_id: str, seq_type: str, mask_feature: str = '1'
    ) -> str:
        """Load Transcript sequence from Ensembl

        Args:
            transcript_id (str): Transcript ID
            seq_type (str): Choice between: genomic,cds,cdna,protein
            mask_feature (str): Masking options for genomic sequence. "0" for no masking, "1" for soft-masking introns & UTRs

        Returns:
            str: Sequence of the transcript
        """
        base_id, _ = self.split_transcript_id(transcript_id)
        ext = f'/sequence/id/{base_id}?type={seq_type};mask_feature={mask_feature}'
        logger.debug(
            f'Fetching sequence for transcript_id: {transcript_id}, seq_type: {seq_type} from Ensembl API.'
        )
        logger.debug(f'Request URL: {self.server + ext}')
        r = self.session.get(
            self.server + ext, headers={'Content-Type': 'text/plain'}, timeout=10
        )
        r.raise_for_status()
        return r.text

    def get_genomic_sequence(
        self,
        chromosome: str,
        start: int,
        end: int,
        strand: str = '1',
        mask_feature: str = '1',
    ) -> str:
        """Load genomic sequence from Ensembl

        Args:
            chromosome (str): Chromosome number (e.g., "1", "X", "Y")
            start (int): Start position (1-based)
            end (int): End position (1-based)
            strand (str): Strand information ("1" for forward, "-1" for reverse)
            mask_feature (str): Masking options for genomic sequence. "0" for no masking, "1" for soft-masking introns & UTRs
        Returns:
            str: Genomic sequence for the specified region
        """
        ext = f'/sequence/region/human/{chromosome}:{start}..{end}:{strand}?mask_feature={mask_feature}'
        logger.debug(
            f'Fetching genomic sequence for region: {chromosome}:{start}-{end}:{strand} from Ensembl API.'
        )
        logger.debug(f'Request URL: {self.server + ext}')
        r = self.session.get(
            self.server + ext, headers={'Content-Type': 'text/plain'}, timeout=10
        )
        r.raise_for_status()
        return r.text

    def map_coordinates(self, transcript_id: str, start: int, end: int, reference: str):
        base_id, _ = self.split_transcript_id(transcript_id)
        # Transform coordinates to 1-based
        start = int(start) + 1
        end = int(end) + 1

        ext = f'/map/{reference}/{base_id}/{start}..{end}'
        r = self.session.get(
            self.server + ext, headers={'Content-Type': 'application/json'}, timeout=10
        )
        r.raise_for_status()
        return r.json()

    def get_transcripts_for_gene(self, gene_id: str) -> list[str]:
        """Get transcript IDs for all transcripts for a given gene (gene ID)

        Args:
            gene_id (str): Gene ID (e.g., "ENSG00000139618")

        Raises:
            ValueError: If the provided gene_id is not an Ensembl Gene ID (does not contain "ENSG")

        Returns:
            list[str]: List of transcript IDs associated with the given gene ID
        """
        if 'ENSG' not in gene_id:
            raise ValueError(
                f'Provided gene_id is not an Ensembl Gene ID. Received: {gene_id}'
            )

        ext = f'/lookup/id/{gene_id}?expand=1'
        r = self.session.get(
            self.server + ext, headers={'Content-Type': 'application/json'}, timeout=10
        )
        r.raise_for_status()
        transcript_ids = [t['id'] for t in r.json().get('Transcript', [])]
        return transcript_ids

    def get_gene_symbol_for_geneID(self, gene_id: str) -> str:
        if 'ENSG' not in gene_id:
            raise ValueError(
                f'Provided gene_id is not an Ensembl Gene ID. Received: {gene_id}'
            )

        ext = f'/lookup/id/{gene_id}'
        r = self.session.get(
            self.server + ext, headers={'Content-Type': 'application/json'}, timeout=10
        )
        r.raise_for_status()
        return r.json().get('display_name', '')

    def get_gene_symbol_for_transcriptID(self, transcript_id: str) -> tuple[str, str]:
        if 'ENST' not in transcript_id:
            raise InvalidTranscriptIdError(
                f'Provided transcript_id is not an Ensembl Transcript ID. Received: {transcript_id}'
            )
        base_id, requested_version = self.split_transcript_id(transcript_id)

        ext = f'/lookup/id/{base_id}'

        r = self.session.get(
            self.server + ext, headers={'Content-Type': 'application/json'}, timeout=10
        )
        r.raise_for_status()

        data = r.json()
        returned_version = data.get('version')
        full_transcript_id = (
            f'{base_id}.{returned_version}' if returned_version else base_id
        )

        if requested_version is not None and returned_version != requested_version:
            raise InvalidTranscriptVersionError(
                f'Transcript version mismatch: requested {base_id}.{requested_version}, '
                f'but Ensembl returned version {base_id}.{returned_version}.'
            )

        # return r.json().get("display_name", "")
        gene_ID = r.json().get('Parent', '')
        gene_symbol = self.get_gene_symbol_for_geneID(gene_ID)
        return gene_symbol, full_transcript_id

    def get_overlapped_geneIDs_for_region(self, chromosome: str, start: int, end: int):
        ext = f'/overlap/region/human/{chromosome}:{start}-{end}?feature=gene'
        r = self.session.get(
            self.server + ext, headers={'Content-Type': 'application/json'}, timeout=10
        )
        r.raise_for_status()
        gene_ids = [gene['id'] for gene in r.json()]
        return gene_ids

    def get_overlapped_geneSymbols_for_region(
        self, chromosome: str, start: int, end: int
    ):
        gene_ids = self.get_overlapped_geneIDs_for_region(chromosome, start, end)

        ext = '/lookup/id'
        r = self.session.post(
            self.server + ext,
            headers={'Content-Type': 'application/json'},
            json={'ids': gene_ids},
            timeout=10,
        )
        r.raise_for_status()
        gene_symbols = [
            gene_info.get('display_name', '') for gene_info in r.json().values()
        ]
        # Drop empty gene symbols
        gene_symbols = [symbol for symbol in gene_symbols if symbol]
        return gene_symbols

    def get_overlapped_genes_details_for_region(
        self, chromosome: str, start: int, end: int
    ):
        gene_ids = self.get_overlapped_geneIDs_for_region(chromosome, start, end)
        gene_symbols = {}
        for gene_id in gene_ids:
            gene_symbols[gene_id] = self.get_gene_symbol_for_geneID(gene_id)
        return gene_symbols
