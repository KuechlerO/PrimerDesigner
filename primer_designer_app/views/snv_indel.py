from django.shortcuts import render, redirect
from django.http import HttpResponse

import logging

from primer_designer_app.models import DesignResultsSummary

from primer_designer_app.utils.helpers import (
    html_visualize_sequence,
    create_hgvs_notation,
    vcf_hits_json_for_display,
    snp_hits_json_for_display,
)
from primer_designer_app.utils.doc_utils import create_primer_report
from primer_designer_app.utils.insilico_analysis import insilico_reference_description
from primer_designer_app.utils.design_validation import validate_primer_search_results
from primer_designer_app.utils.variant_info import SequenceVariantInfo
from primer_designer_app.views.view_utils import (
    _get_post,
    build_form_data_from_request,
    build_primer_settings,
    handle_transcript_input,
    handle_genomic_snv,
    handle_genomic_indel,
    handle_sequence_input,
)

logger = logging.getLogger(__name__)


def index(request):
    """
    View to render the idex page for the snv-indels page
    """
    return render(
        request,
        "primer_designer_app/snv_indel_index.html",
        {"form_data": build_form_data_from_request(request)},
    )


def primers_overview(request, uuid=None):
    """Generates the output page with primer overview and highlighted sequence snippet."""
    logger.debug(
        f"Received request for primers overview with UUID: {uuid} and POST data: {request.POST}"
    )

    # Fetch design results object using the provided UUID
    if uuid is None:
        # Build primer settings object from POST data
        primer_settings_obj = build_primer_settings(request)

        try:
            # --- Handle transcript ID input first ---
            if _get_post(request, "Transcript-ID", None):
                new_uuid = handle_transcript_input(request, primer_settings_obj)

            # --- Handle genomic input next ---
            # SNV (genomic) path
            elif _get_post(request, "genom_pos", None):
                new_uuid = handle_genomic_snv(request, primer_settings_obj)

            elif all(
                _get_post(request, input, None)
                for input in ["IndelChrom", "IndelStart", "IndelEnd", "IndelIns"]
            ):
                new_uuid = handle_genomic_indel(request, primer_settings_obj)

            # --- Handle sequence input last ---
            elif _get_post(request, "sequence", None):
                new_uuid = handle_sequence_input(request, primer_settings_obj)
            else:
                return HttpResponse(
                    "Invalid input: No recognizable input field found.", status=400
                )
        except ValueError as exc:
            logger.warning("Invalid primer design input: %s", exc)
            return HttpResponse(str(exc), status=400)

        designResults_obj = DesignResultsSummary.objects.get(id=new_uuid)

        logger.debug(f"Created new DesignResultsSummary object with UUID: {new_uuid}")
        logger.debug(
            f"DesignResultsSummary variant info: {designResults_obj.get_variant_info()}"
        )

    else:
        designResults_obj = DesignResultsSummary.objects.get(id=uuid)
        logger.debug(f"DesignResultsSummary data2: {designResults_obj}")

    prim_search_results = designResults_obj.get_primer_search_results()
    validate_primer_search_results(prim_search_results)
    var_info = designResults_obj.get_variant_info()

    highlighted_seq_snippet, display_offset, display_length, display_chunks = (
        html_visualize_sequence(
            designResults_obj.primer_settings,
            var_info,
            prim_search_results.primer_pairs[0],
            all_primer_pairs=prim_search_results.primer_pairs,
        )
    )

    primerF_sequences, primerR_sequences = (
        [primer_pair.left_seq for primer_pair in prim_search_results.primer_pairs],
        [primer_pair.right_seq for primer_pair in prim_search_results.primer_pairs],
    )

    hgvs_info = create_hgvs_notation(var_info)

    logger.debug(f"primerF_sequences: {primerF_sequences}")

    # Print amplicons per primer pair for debugging
    for i, primer_pair in enumerate(prim_search_results.primer_pairs):
        logger.debug(f"Primer pair {i}: Amplicons per primer: {primer_pair.amplicons}")

    snp_analysis = designResults_obj.snp_analysis_data or {}
    snp_hits_json = "[]"
    if snp_analysis.get("hits"):
        snp_hits_json = snp_hits_json_for_display(
            snp_analysis["hits"], display_offset, display_length
        )

    vcf_applied = (designResults_obj.variant_info_data or {}).get(
        "vcf_applied_variants"
    ) or []
    vcf_hits_json = vcf_hits_json_for_display(
        vcf_applied, display_offset, display_length
    )

    is_sequence_only = isinstance(var_info, SequenceVariantInfo) and not getattr(
        var_info, "genomic_pos", None
    )
    amplicon_requested = bool(
        getattr(designResults_obj.primer_settings, "do_insilico_pcr", False)
    )
    show_amplicon_sequence_warning = is_sequence_only and amplicon_requested
    show_snp_sequence_warning = (
        bool(snp_analysis.get("enabled"))
        and is_sequence_only
        and "requires genomic or transcript coordinates"
        in str(snp_analysis.get("message", ""))
    )

    return render(
        request,
        "primer_designer_app/snv_indel_results.html",
        {
            "highlighted_sequence": highlighted_seq_snippet,
            "display_chunks": display_chunks,
            "result_sum_obj": designResults_obj,
            "hgvs_info": hgvs_info,
            "snp_analysis": snp_analysis,
            "snp_hits_json": snp_hits_json,
            "vcf_hits_json": vcf_hits_json,
            "sequence_display_offset": display_offset,
            "show_amplicon_sequence_warning": show_amplicon_sequence_warning,
            "show_snp_sequence_warning": show_snp_sequence_warning,
            "insilico_reference_note": insilico_reference_description(
                designResults_obj.primer_settings
            ),
        },
    )


def primer_details(request, uuid):
    """
    View to display the final outut after the user has decided on a primer pair.
    """

    def generate_output(request, uuid, selected_primer_index: int):
        # Fetch saved object (infos and results) with UUID
        retrieved_result = DesignResultsSummary.objects.get(id=uuid)
        logger.debug(f"Retrieved DesignResultsSummary object with UUID: {uuid}")
        logger.debug(
            f"Retrieved primer search results: {retrieved_result.primer_search_results}"
        )
        logger.debug(f"Retrieved variant info: {retrieved_result.get_variant_info()}")

        primer_pairs = retrieved_result.get_primer_search_results().primer_pairs
        selected_primer = primer_pairs[selected_primer_index - 1]

        var_info_detail = retrieved_result.get_variant_info()
        full_highlighted_seq, display_offset, display_length, display_chunks = (
            html_visualize_sequence(
                retrieved_result.primer_settings,
                var_info_detail,
                selected_primer,
            )
        )

        genome_version = retrieved_result.get_variant_info().ref_genome
        hgvs_notation = create_hgvs_notation(retrieved_result.get_variant_info())

        snp_analysis = retrieved_result.snp_analysis_data or {}
        snp_hits_json = "[]"
        if snp_analysis.get("hits"):
            snp_hits_json = snp_hits_json_for_display(
                snp_analysis["hits"], display_offset, display_length
            )

        vcf_applied = (retrieved_result.variant_info_data or {}).get(
            "vcf_applied_variants"
        ) or []
        vcf_hits_json = vcf_hits_json_for_display(
            vcf_applied, display_offset, display_length
        )

        return render(
            request,
            "primer_designer_app/snv_indel_details.html",
            {
                "highlighted_sequence": full_highlighted_seq,
                "display_chunks": display_chunks,
                "primer_pair": selected_primer,
                "uuid": retrieved_result.id,
                "genome_version": genome_version,
                "hgvs_info": hgvs_notation,
                "selected_primer_index": selected_primer_index,
                "snp_analysis": snp_analysis,
                "snp_hits_json": snp_hits_json,
                "vcf_hits_json": vcf_hits_json,
                "sequence_display_offset": display_offset,
                "insilico_reference_note": insilico_reference_description(
                    retrieved_result.primer_settings
                ),
            },
        )

    # process only after POST rerquest
    if request.method == "POST":
        selected_primer_index = int(request.POST.get("selected-primer", None))
        # fail guard
        if selected_primer_index is None:
            return HttpResponse("Kein Primer ausgewählt.", status=400)

        return generate_output(request, uuid, selected_primer_index)
    # else redirect to the overview page
    return redirect("primer_designer_app:snv_indel_primers_overview", uuid=uuid)


def generate_report(request, uuid, selected_primer_index):
    """
    View to generate the final report after the user has decided on a primer pair.
    """
    designResults_obj = DesignResultsSummary.objects.get(id=uuid)
    doc_buffer = create_primer_report(designResults_obj, selected_primer_index)
    response = HttpResponse(
        doc_buffer,
        content_type="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
    )
    response["Content-Disposition"] = (
        f"attachment; filename=report_{uuid}_{selected_primer_index}.docx"
    )
    return response
