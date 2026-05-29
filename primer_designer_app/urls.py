from django.urls import path

from .views import snv_indel
from .views import structural_variant
from .views import documentation_view
from .views import allele_specific

app_name = "primer_designer_app"

urlpatterns = [
    path("", snv_indel.index, name="index"),
    # ---- SNVs and Indels ----
    # Index
    path("snv-indel/", snv_indel.index, name="snv_indel_index"),
    # Overview
    path(
        "snv-indel/primers-overview/",
        snv_indel.primers_overview,
        name="snv_indel_primers_overview",
    ),
    path(
        "snv-indel/primers-overview/<uuid:uuid>/",
        snv_indel.primers_overview,
        name="snv_indel_primers_overview_with_uuid",
    ),
    # Details view
    path(
        "snv-indel/primer-details/<uuid:uuid>/",
        snv_indel.primer_details,
        name="snv_indel_primer_details",
    ),
    # Report creation
    path(
        "snv-indel/generate-report/<uuid:uuid>/<int:selected_primer_index>/",
        snv_indel.generate_report,
        name="snv_indel_generate_report",
    ),
    # ---- Structural Variants ----
    path(
        "structural-variant/",
        structural_variant.index,
        name="structural_variants_index",
    ),
    path(
        "structural-variant/generate-report/<uuid:uuid>/",
        structural_variant.generate_report,
        name="structural_variant_generate_report",
    ),
    # ---- Allele-specific PCR (AS-PCR) ----
    path("allele-specific/", allele_specific.index, name="allele_specific_index"),
    path(
        "allele-specific/primers-overview/",
        allele_specific.primers_overview,
        name="allele_specific_primers_overview",
    ),
    path(
        "allele-specific/primers-overview/<uuid:uuid>/",
        allele_specific.primers_overview,
        name="allele_specific_primers_overview_with_uuid",
    ),
    path(
        "allele-specific/primer-details/<uuid:uuid>/",
        allele_specific.primer_details,
        name="allele_specific_primer_details",
    ),
    path(
        "allele-specific/generate-report/<uuid:uuid>/<int:selected_primer_index>/",
        allele_specific.generate_report,
        name="allele_specific_generate_report",
    ),
    # Documentation
    path(
        "documentation/",
        documentation_view.documentation,
        name="documentation",
    ),
]
