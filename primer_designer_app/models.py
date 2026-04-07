# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey and OneToOneField has `on_delete` set to the desired behavior
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from enum import Enum
from django.db import models
from django.db import connection

from primer_designer_app.utils.variant_info import (
    VariantInfo,
    TranscriptVariantInfo,
    GenomicVariantInfo,
    SequenceVariantInfo,
    IndelType,
    ReferenceType,
)
from primer_designer_app.utils.primer_utils import PrimerSearchResults

import uuid
import logging

LOGGER = logging.getLogger(__name__)


class PrimerSettingsModel(models.Model):
    # Bases on each side of the variant included in Primer3 SEQUENCE_TARGET (was PCR=50, qPCR=30).
    target_padding = models.IntegerField(default=50)
    tm = models.IntegerField()
    gc = models.IntegerField()
    max_poly_x = models.IntegerField()
    productsize_range = models.JSONField()  # Store as a list in JSON format
    reference_genome = models.CharField(max_length=50)
    primer_size = models.IntegerField(null=True, blank=True)
    target = models.JSONField(null=True, blank=True)  # Store as a list [start, length]
    # Context for insilico analysis: Either "transcriptomic" or "genomic"
    context = models.CharField(
        max_length=20, null=False, blank=False, default='genomic'
    )
    # Run Dicey in-silico PCR / amplicon search (optional; default off)
    do_insilico_pcr = models.BooleanField(default=False)
    # Optional Primer3 global_args overrides (custom mode); merged in primer_utils
    primer3_overrides = models.JSONField(default=dict, blank=True)

    def set_target(self, rel_pos):
        LOGGER.debug(
            'Reference genome: %s, target_padding: %s, relative position: %s',
            self.reference_genome,
            self.target_padding,
            rel_pos,
        )
        """Set Primer3 SEQUENCE_TARGET from variant interval and per-side padding (bp)."""
        offset = int(self.target_padding)
        if offset < 1 or offset > 500:
            raise ValueError(f'Invalid target_padding: {offset}')

        nr_deleted_bases = rel_pos[1] - rel_pos[0]
        self.target = [rel_pos[0] - offset, nr_deleted_bases + 2 * offset]
        self.save()

    def set_context(self, context):
        """Set genomic vs transcriptomic context for in-silico (Dicey) analysis."""
        if context in ['transcriptomic', 'genomic']:
            self.context = context
            self.save()
        else:
            raise ValueError(
                f"Invalid context: {context}. Must be 'transcriptomic' or 'genomic'."
            )


class DesignResultsSummary(models.Model):
    """Model to store primer design parameters and results in a JSONField for later retrieval."""

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    created_at = models.DateTimeField(auto_now_add=True)

    primer_settings = models.OneToOneField(
        'PrimerSettingsModel',  # Replace with the actual model name for PrimerSettings
        on_delete=models.CASCADE,
        null=True,
        related_name='param_forwarding',
    )

    variant_info_data = models.JSONField(
        null=True, blank=True
    )  # Store serialized VariantInfo
    primer_search_results = models.JSONField(
        null=True, blank=True
    )  # Store primer results as JSON

    def save_primer_results(self, prim_search_res, prim_settings, var_info):
        """
        Save primer results and related information into the DesignResultsSummary model
        """

        # Save references to the related objects
        self.primer_settings = prim_settings
        self.variant_info_data = (
            var_info.__dict__
        )  # Serialize VariantInfo to a dictionary
        self.primer_search_results = (
            prim_search_res.__dict__
        )  # Serialize PrimerSearchResults to a dictionary

        # Serialize VariantInfo to a dictionary and convert non-serializable fields
        variant_info_dict = var_info.__dict__.copy()
        if isinstance(variant_info_dict.get('indel_type'), Enum):
            variant_info_dict['indel_type'] = variant_info_dict[
                'indel_type'
            ].value  # Convert Enum to string
        if isinstance(variant_info_dict.get('reference_type'), Enum):
            variant_info_dict['reference_type'] = variant_info_dict[
                'reference_type'
            ].value  # Convert Enum to string
        self.variant_info_data = variant_info_dict

        # Serialize PrimerSearchResults to a dictionary
        primer_results_dict = prim_search_res.__dict__.copy()
        if 'primer_pairs' in primer_results_dict:
            # Convert PrimerPairResult objects to dictionaries
            primer_results_dict['primer_pairs'] = [
                pair.to_dict() for pair in primer_results_dict['primer_pairs']
            ]
        self.primer_search_results = primer_results_dict

        # Save the DesignResultsSummary instance
        self.save()

        LOGGER.debug(f"Saved DesignResultsSummary with ID: {self.id}")
        LOGGER.debug(f"VariantInfo saved: {self.variant_info_data}")
        LOGGER.debug(f"PrimerSearchResults saved: {self.primer_search_results}")
        # If the current DB backend is not PostgreSQL, skip creating Postgres-only
        # normalized rows (ArrayField etc.) to avoid binding errors (e.g. SQLite).
        if connection.vendor != 'postgresql':
            # normalized models require Postgres (ArrayField); nothing more to do
            return

    def get_variant_info(self):
        """Deserialize the stored variant info from JSON back into a VariantInfo object."""
        if self.variant_info_data:
            # Convert IndelType from string back to Enum if necessary
            if 'indel_type' in self.variant_info_data:
                self.variant_info_data['indel_type'] = IndelType(
                    self.variant_info_data['indel_type']
                )
            if 'reference_type' in self.variant_info_data:
                self.variant_info_data['reference_type'] = ReferenceType(
                    self.variant_info_data['reference_type']
                )

            LOGGER.debug(
                f"Deserializing VariantInfo from data: {self.variant_info_data}"
            )
            # Test if TranscriptVariantInfo or else VariantInfo
            if 'transcript_id' in self.variant_info_data:
                return TranscriptVariantInfo(**self.variant_info_data)
            elif (
                'genomic_pos' in self.variant_info_data
                and self.variant_info_data['genomic_pos']
            ):
                LOGGER.debug(
                    f"Deserializing GenomicVariantInfo from data: {self.variant_info_data}"
                )
                LOGGER.debug(
                    f"Genomic position: {self.variant_info_data.get('genomic_pos')}"
                )
                return GenomicVariantInfo(**self.variant_info_data)
            else:
                # First create a VariantInfo object
                tmp_varInfo_obj = VariantInfo(**self.variant_info_data)
                input_seq = tmp_varInfo_obj.get_seq('input')
                return SequenceVariantInfo(
                    input_seq=input_seq, **self.variant_info_data
                )
        return None

    def get_primer_search_results(self):
        """Deserialize the stored primer search results from JSON back into a PrimerSearchResults object."""
        if self.primer_search_results:
            return PrimerSearchResults.from_dict(self.primer_search_results)
        return None
