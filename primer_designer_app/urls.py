from django.urls import path

from .views import snv_indel
from .views import structural_variant


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
    path(
        "structural-variant/",
        structural_variant.index,
        name="structural_variants_index",
    ),
]
