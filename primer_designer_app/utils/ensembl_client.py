import logging
import json
import time

import requests
from requests.adapters import HTTPAdapter, Retry

from primer_designer_app.exceptions import (
    InvalidTranscriptIdError,
    InvalidTranscriptVersionError,
)

server37 = "https://grch37.rest.ensembl.org"
server38 = "https://rest.ensembl.org"

logger = logging.getLogger(__name__)

# Ensembl POST /variation allows at most 200 IDs per request.
VARIATION_BATCH_SIZE = 200
GNOMAD_VARIANT_SET = "gnomAD"

# (connect, read) timeouts in seconds
TIMEOUT_LOOKUP = (5, 30)
TIMEOUT_SEQUENCE = (5, 90)
TIMEOUT_OVERLAP = (5, 30)
TIMEOUT_VARIATION = (5, 60)
# rest.ensembl.org often returns sticky short-lived 500s; wait and retry.
_SEQUENCE_RETRY_DELAYS_SEC = (0.0, 1.5, 3.0, 6.0, 12.0)

# #region agent log
_DEBUG_LOG_PATH = "/Users/oliverkuchler/Programming/git_projects/PrimerDesigner/.cursor/debug-e4e356.log"


def _agent_log(hypothesis_id: str, location: str, message: str, data: dict | None = None):
    try:
        payload = {
            "sessionId": "e4e356",
            "hypothesisId": hypothesis_id,
            "location": location,
            "message": message,
            "data": data or {},
            "timestamp": int(time.time() * 1000),
        }
        with open(_DEBUG_LOG_PATH, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(payload) + "\n")
    except Exception:
        pass


# #endregion


class EnsemblClient:
    def __init__(self, ref_genome: str = "GRCh38"):
        if ref_genome == "GRCh38":
            self.server = server38
        elif ref_genome == "GRCh37":
            self.server = server37
        else:
            raise ValueError(
                f"Unsupported reference genome: {ref_genome}. Supported values are 'GRCh38' and 'GRCh37'."
            )

        self.session = requests.Session()
        # Retry rate-limits and gateway errors, but not generic 500s — rest.ensembl.org
        # often returns sticky 500s where retries only amplify Max retries exceeded.
        retries = Retry(
            total=2,
            backoff_factor=0.8,
            status_forcelist=[429, 502, 503, 504],
            allowed_methods=frozenset(["GET", "POST"]),
        )
        adapter = HTTPAdapter(max_retries=retries)
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)

    def split_transcript_id(self, transcript_id: str):
        transcript_id = transcript_id.strip()

        if "." in transcript_id:
            base_id, version = transcript_id.split(".", 1)
            return base_id, int(version)

        return transcript_id, None

    def get_transcript_sequence(
        self, transcript_id: str, seq_type: str, mask_feature: str | None = None
    ) -> str:
        """Load Transcript sequence from Ensembl

        Args:
            transcript_id (str): Transcript ID
            seq_type (str): Choice between: genomic,cds,cdna,protein
            mask_feature (str | None): Only for genomic-style masking. Omit for
                cdna/cds — ``mask_feature=0/1`` currently 500s on rest.ensembl.org.

        Returns:
            str: Sequence of the transcript
        """
        base_id, _ = self.split_transcript_id(transcript_id)
        ext = f"/sequence/id/{base_id}?type={seq_type}"
        if mask_feature is not None:
            ext += f";mask_feature={mask_feature}"
        url = self.server + ext
        logger.debug(
            f"Fetching sequence for transcript_id: {transcript_id}, seq_type: {seq_type} from Ensembl API."
        )
        logger.debug(f"Request URL: {url}")
        # #region agent log
        _agent_log(
            "A",
            "ensembl_client.py:get_transcript_sequence:before",
            "Starting transcript sequence GET",
            {
                "server": self.server,
                "url": url,
                "seq_type": seq_type,
                "mask_feature": mask_feature,
                "transcript_id": transcript_id,
            },
        )
        # #endregion
        t0 = time.time()
        try:
            seq = self._get_sequence_text(url)
            elapsed_ms = int((time.time() - t0) * 1000)
            # #region agent log
            _agent_log(
                "A",
                "ensembl_client.py:get_transcript_sequence:after",
                "Transcript sequence GET completed",
                {
                    "server": self.server,
                    "url": url,
                    "status_code": 200,
                    "elapsed_ms": elapsed_ms,
                    "body_len": len(seq or ""),
                    "body_prefix": (seq or "")[:40],
                },
            )
            # #endregion
            return seq
        except Exception as exc:
            elapsed_ms = int((time.time() - t0) * 1000)
            # #region agent log
            _agent_log(
                "A",
                "ensembl_client.py:get_transcript_sequence:error",
                "Transcript sequence GET failed",
                {
                    "server": self.server,
                    "url": url,
                    "elapsed_ms": elapsed_ms,
                    "error_type": type(exc).__name__,
                    "error": str(exc)[:500],
                },
            )
            # #endregion
            raise

    def _get_sequence_text(self, url: str) -> str:
        """Fetch a sequence endpoint with plain/JSON attempts and delayed 500 retries.

        rest.ensembl.org frequently returns immediate HTML 500 pages under load, then
        succeeds after a few seconds (successful reads can take 20–40s).
        """
        last_error: Exception | None = None
        for attempt, delay in enumerate(_SEQUENCE_RETRY_DELAYS_SEC):
            if delay:
                time.sleep(delay)
            for headers, mode in (
                ({"Content-Type": "text/plain", "Accept": "text/plain"}, "plain"),
                (
                    {
                        "Content-Type": "application/json",
                        "Accept": "application/json",
                    },
                    "json",
                ),
            ):
                try:
                    r = self.session.get(
                        url, headers=headers, timeout=TIMEOUT_SEQUENCE
                    )
                    if r.status_code == 200:
                        if mode == "plain":
                            text = r.text or ""
                            if text and not text.lstrip().startswith("<!"):
                                # #region agent log
                                _agent_log(
                                    "A",
                                    "ensembl_client.py:_get_sequence_text:success",
                                    "Sequence fetch succeeded",
                                    {
                                        "url": url,
                                        "mode": mode,
                                        "attempt": attempt,
                                        "body_len": len(text),
                                    },
                                )
                                # #endregion
                                return text
                        else:
                            payload = r.json()
                            if isinstance(payload, dict) and payload.get("seq"):
                                # #region agent log
                                _agent_log(
                                    "A",
                                    "ensembl_client.py:_get_sequence_text:success",
                                    "Sequence fetch succeeded",
                                    {
                                        "url": url,
                                        "mode": mode,
                                        "attempt": attempt,
                                        "body_len": len(payload["seq"]),
                                    },
                                )
                                # #endregion
                                return payload["seq"]

                    # #region agent log
                    _agent_log(
                        "C",
                        "ensembl_client.py:_get_sequence_text:attempt",
                        "Sequence attempt failed",
                        {
                            "url": url,
                            "mode": mode,
                            "attempt": attempt,
                            "status_code": r.status_code,
                            "body_prefix": (r.text or "")[:50],
                        },
                    )
                    # #endregion
                    if r.status_code >= 500:
                        last_error = requests.HTTPError(
                            f"{r.status_code} Server Error for url: {url}",
                            response=r,
                        )
                        continue
                    r.raise_for_status()
                except Exception as exc:
                    last_error = exc
                    # #region agent log
                    _agent_log(
                        "C",
                        "ensembl_client.py:_get_sequence_text:attempt_exc",
                        "Sequence attempt raised",
                        {
                            "url": url,
                            "mode": mode,
                            "attempt": attempt,
                            "error_type": type(exc).__name__,
                            "error": str(exc)[:300],
                        },
                    )
                    # #endregion
                    continue

        if last_error is not None:
            raise last_error
        raise ValueError(f"Unexpected Ensembl sequence response for {url}")

    def get_genomic_sequence(
        self,
        chromosome: str,
        start: int,
        end: int,
        strand: str = "1",
        mask_feature: str = "1",
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
        ext = f"/sequence/region/human/{chromosome}:{start}..{end}:{strand}?mask_feature={mask_feature}"
        logger.debug(
            f"Fetching genomic sequence for region: {chromosome}:{start}-{end}:{strand} from Ensembl API."
        )
        logger.debug(f"Request URL: {self.server + ext}")
        r = self.session.get(
            self.server + ext,
            headers={"Content-Type": "text/plain"},
            timeout=TIMEOUT_SEQUENCE,
        )
        r.raise_for_status()
        return r.text

    def map_coordinates(self, transcript_id: str, start: int, end: int, reference: str):
        base_id, _ = self.split_transcript_id(transcript_id)
        # Transform coordinates to 1-based
        start = int(start) + 1
        end = int(end) + 1

        ext = f"/map/{reference}/{base_id}/{start}..{end}"
        r = self.session.get(
            self.server + ext,
            headers={"Content-Type": "application/json"},
            timeout=TIMEOUT_LOOKUP,
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
        if "ENSG" not in gene_id:
            raise ValueError(
                f"Provided gene_id is not an Ensembl Gene ID. Received: {gene_id}"
            )

        ext = f"/lookup/id/{gene_id}?expand=1"
        r = self.session.get(
            self.server + ext,
            headers={"Content-Type": "application/json"},
            timeout=TIMEOUT_LOOKUP,
        )
        r.raise_for_status()
        transcript_ids = [t["id"] for t in r.json().get("Transcript", [])]
        return transcript_ids

    def get_gene_symbol_for_geneID(self, gene_id: str) -> str:
        if "ENSG" not in gene_id:
            raise ValueError(
                f"Provided gene_id is not an Ensembl Gene ID. Received: {gene_id}"
            )

        ext = f"/lookup/id/{gene_id}"
        r = self.session.get(
            self.server + ext,
            headers={"Content-Type": "application/json"},
            timeout=TIMEOUT_LOOKUP,
        )
        r.raise_for_status()
        return r.json().get("display_name", "")

    @staticmethod
    def gene_symbol_from_transcript_display_name(display_name: str) -> str:
        """Map Ensembl transcript display names like ``PLPP5-203`` to ``PLPP5``."""
        if not display_name:
            return ""
        return str(display_name).split("-", 1)[0].strip()

    def get_current_gene_symbol(self, gene_id: str) -> str:
        """Return the gene display_name from current Ensembl (GRCh38), if available.

        GRCh37 archive still uses retired symbols (e.g. PPAPDC1B vs PLPP5).
        """
        if "ENSG" not in gene_id:
            return ""
        ext = f"/lookup/id/{gene_id}"
        try:
            r = self.session.get(
                server38 + ext,
                headers={"Content-Type": "application/json"},
                timeout=TIMEOUT_LOOKUP,
            )
            r.raise_for_status()
            return r.json().get("display_name", "") or ""
        except Exception:
            logger.debug(
                "Could not resolve current gene symbol for %s via %s",
                gene_id,
                server38,
                exc_info=True,
            )
            return ""

    def get_gene_symbol_for_transcriptID(
        self, transcript_id: str
    ) -> tuple[str, str, str]:
        """Return ``(gene_symbol, gene_ID, full_transcript_id)``.

        Uses the transcript lookup first (including ``display_name``). A follow-up
        gene lookup is best-effort because ``rest.ensembl.org`` often 500s on
        ``/lookup/id/ENSG...`` even when the transcript lookup succeeds.
        """
        if "ENST" not in transcript_id:
            raise InvalidTranscriptIdError(
                f"Provided transcript_id is not an Ensembl Transcript ID. Received: {transcript_id}"
            )
        base_id, requested_version = self.split_transcript_id(transcript_id)

        ext = f"/lookup/id/{base_id}"

        t0 = time.time()
        try:
            r = self.session.get(
                self.server + ext,
                headers={"Content-Type": "application/json"},
                timeout=TIMEOUT_LOOKUP,
            )
            # #region agent log
            _agent_log(
                "B",
                "ensembl_client.py:get_gene_symbol_for_transcriptID",
                "Transcript lookup response",
                {
                    "server": self.server,
                    "url": self.server + ext,
                    "status_code": r.status_code,
                    "elapsed_ms": int((time.time() - t0) * 1000),
                },
            )
            # #endregion
            r.raise_for_status()
        except Exception as exc:
            # #region agent log
            _agent_log(
                "B",
                "ensembl_client.py:get_gene_symbol_for_transcriptID:error",
                "Transcript lookup failed",
                {
                    "server": self.server,
                    "url": self.server + ext,
                    "elapsed_ms": int((time.time() - t0) * 1000),
                    "error_type": type(exc).__name__,
                    "error": str(exc)[:500],
                },
            )
            # #endregion
            raise

        data = r.json()
        returned_version = data.get("version")
        full_transcript_id = (
            f"{base_id}.{returned_version}" if returned_version else base_id
        )

        if requested_version is not None and returned_version != requested_version:
            raise InvalidTranscriptVersionError(
                f"Transcript version mismatch: requested {base_id}.{requested_version}, "
                f"but Ensembl returned version {base_id}.{returned_version}."
            )

        gene_ID = data.get("Parent", "") or ""
        gene_symbol = self.gene_symbol_from_transcript_display_name(
            data.get("display_name", "")
        )
        if gene_ID:
            try:
                looked_up = self.get_gene_symbol_for_geneID(gene_ID)
                if looked_up:
                    gene_symbol = looked_up
            except Exception as exc:
                # #region agent log
                _agent_log(
                    "E",
                    "ensembl_client.py:get_gene_symbol_for_transcriptID:gene_lookup",
                    "Gene lookup failed; keeping transcript-derived symbol",
                    {
                        "server": self.server,
                        "gene_id": gene_ID,
                        "fallback_symbol": gene_symbol,
                        "error_type": type(exc).__name__,
                        "error": str(exc)[:300],
                    },
                )
                # #endregion
                logger.warning(
                    "Ensembl gene lookup failed for %s (%s); using transcript symbol %s",
                    gene_ID,
                    exc,
                    gene_symbol,
                )
            # Prefer modern HGNC/Ensembl symbol when designing on GRCh37 archive.
            if self.server != server38:
                current = self.get_current_gene_symbol(gene_ID)
                if current:
                    gene_symbol = current
        return gene_symbol, gene_ID, full_transcript_id

    def get_overlapping_variations_for_region(
        self,
        chromosome: str,
        start: int,
        end: int,
        *,
        variant_set: str | None = None,
    ) -> list:
        """Return variation features overlapping a 1-based inclusive genomic region."""
        ext = f"/overlap/region/human/{chromosome}:{start}-{end}?feature=variation"
        if variant_set:
            ext += f";variant_set={variant_set}"
        logger.debug(
            "Fetching variation overlap for %s:%s-%s (variant_set=%s) from Ensembl API.",
            chromosome,
            start,
            end,
            variant_set,
        )
        r = self.session.get(
            self.server + ext,
            headers={"Content-Type": "application/json"},
            timeout=TIMEOUT_OVERLAP,
        )
        r.raise_for_status()
        return r.json()

    def get_variation_details_batch(self, variant_ids: list[str]) -> dict[str, dict]:
        """Fetch population frequencies for rsIDs via POST /variation (chunks of 200)."""
        if not variant_ids:
            return {}

        unique_ids = list(dict.fromkeys(v for v in variant_ids if v))
        results: dict[str, dict] = {}
        ext = "/variation/human?pops=1"

        for offset in range(0, len(unique_ids), VARIATION_BATCH_SIZE):
            chunk = unique_ids[offset : offset + VARIATION_BATCH_SIZE]
            logger.debug(
                "Fetching variation frequencies for %s IDs (batch %s-%s).",
                len(chunk),
                offset,
                offset + len(chunk),
            )
            r = self.session.post(
                self.server + ext,
                headers={"Content-Type": "application/json"},
                json={"ids": chunk},
                timeout=TIMEOUT_VARIATION,
            )
            r.raise_for_status()
            batch = r.json()
            if isinstance(batch, dict):
                results.update(batch)

        return results

    def get_overlapped_geneIDs_for_region(self, chromosome: str, start: int, end: int):
        ext = f"/overlap/region/human/{chromosome}:{start}-{end}?feature=gene"
        r = self.session.get(
            self.server + ext,
            headers={"Content-Type": "application/json"},
            timeout=TIMEOUT_OVERLAP,
        )
        r.raise_for_status()
        gene_ids = [gene["id"] for gene in r.json()]
        return gene_ids

    def get_overlapped_geneSymbols_for_region(
        self, chromosome: str, start: int, end: int
    ):
        gene_ids = self.get_overlapped_geneIDs_for_region(chromosome, start, end)

        ext = "/lookup/id"
        r = self.session.post(
            self.server + ext,
            headers={"Content-Type": "application/json"},
            json={"ids": gene_ids},
            timeout=TIMEOUT_LOOKUP,
        )
        r.raise_for_status()
        gene_symbols = [
            gene_info.get("display_name", "") for gene_info in r.json().values()
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
