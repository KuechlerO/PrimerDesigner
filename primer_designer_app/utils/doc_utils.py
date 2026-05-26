import io
from datetime import datetime
from typing import Optional

from primer_designer_app.models import PrimerSettingsModel, DesignResultsSummary
from primer_designer_app.utils.variant_info import (
    AllelicVariantInfo,
    TranscriptVariantInfo,
)
from primer_designer_app.utils.helpers import create_hgvs_notation
from primer_designer_app.utils.amplicon_display import (
    amplicon_chrom_label,
    extract_amplicon_summary,
    format_penalty_score,
    truncate_product_seq,
)
from primer_designer_app.utils.insilico_analysis import insilico_reference_description
from primer_designer_app.utils.primer_utils import INSILICO_OK, PrimerPairResult
from primer_designer_app.utils.sv_storage import SV_WINDOW_ORDER

from django.conf import settings
from django.urls import reverse

from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

import logging

logger = logging.getLogger(__name__)


def visualize_sequence_as_docx(
    paragraph,
    prim_settings: PrimerSettingsModel,
    var_info: AllelicVariantInfo,
    selected_primer_pair: PrimerPairResult,
    line_length: int = 60,
):

    seq = var_info.get_seq("input")

    primerF_start, primerF_end = [
        selected_primer_pair.left_relPos_start,
        selected_primer_pair.left_relPos_end,
    ]
    r_primers_offset = (
        max(1, len(var_info.ref_bases)) + 3
    )  # +3 for the brackets in the sequence annotation
    primerR_start, primerR_end = [
        selected_primer_pair.right_relPos_start + r_primers_offset,
        selected_primer_pair.right_relPos_end + r_primers_offset,
    ]
    logger.debug(
        "Primer binding sites (relative to target region): "
        "Forward: %s-%s, Reverse: %s-%s",
        primerF_start,
        primerF_end,
        primerR_start,
        primerR_end,
    )

    start_target = max(prim_settings.target[0], primerF_end + 1)
    end_target = min(
        prim_settings.target[0] + prim_settings.target[1], primerR_start - 1
    )

    # Iterate through the sequence
    counter = 0
    in_mutation_region = False
    for i, char in enumerate(seq):
        run = paragraph.add_run(char)
        run.font.name = "Courier New"  # Schriftart auf 'PT Mono' ändern
        run.font.size = Pt(11)
        if (i >= primerF_start and i <= primerF_end) or (
            i >= primerR_start and i <= primerR_end
        ):
            # Start of forward primer binding site
            # Highlight this region in the document
            run.font.highlight_color = 3
        elif i >= start_target and i < end_target:
            if char == "[":
                in_mutation_region = True
            elif char == "]":
                in_mutation_region = False

            if in_mutation_region or char in ["[", "]"]:
                run.font.highlight_color = 4
            else:
                # Highlight target region
                run.font.highlight_color = 5

        counter += 1
        if counter >= line_length:
            paragraph.add_run(f"  {format(i+1, ',')}\n")
            counter = 0


def add_hyperlink(paragraph, url, text, color="0000FF", underline=True):
    # Create the w:hyperlink tag and add needed values
    part = paragraph.part
    r_id = part.relate_to(
        url,
        "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink",
        is_external=True,
    )
    hyperlink = OxmlElement("w:hyperlink")
    hyperlink.set(qn("r:id"), r_id)

    # Create a w:r element
    new_run = OxmlElement("w:r")

    # Create a w:rPr element
    rPr = OxmlElement("w:rPr")

    # Add color if provided
    if color:
        c = OxmlElement("w:color")
        c.set(qn("w:val"), color)
        rPr.append(c)

    # Underline
    if underline:
        u = OxmlElement("w:u")
        u.set(qn("w:val"), "single")
        rPr.append(u)

    new_run.append(rPr)

    # Create a w:t element and set the text
    text_elem = OxmlElement("w:t")
    text_elem.text = text
    new_run.append(text_elem)

    hyperlink.append(new_run)
    paragraph._p.append(hyperlink)
    return hyperlink


def _set_cell_run(
    cell, text: str, *, bold: bool = False, size_pt: Optional[float] = None
):
    """Replace cell content with a single run (for table cells)."""
    cell.text = ""
    paragraph = cell.paragraphs[0]
    run = paragraph.add_run(str(text))
    if bold:
        run.bold = True
    if size_pt is not None:
        run.font.size = Pt(size_pt)


def add_amplicon_detail_table_to_doc(doc, primer_pair: PrimerPairResult, prim_settings):
    """
    Append a Dicey-style amplicon table when the pair has hits (insilico_status ok).
    """
    if getattr(primer_pair, "insilico_status", None) != INSILICO_OK:
        return
    amplicons = primer_pair.amplicons or []
    if not amplicons:
        return

    doc.add_heading("In-silico amplicons (Dicey)", level=2)
    ref_note = insilico_reference_description(prim_settings)
    if ref_note:
        note_p = doc.add_paragraph(ref_note)
        note_p.paragraph_format.space_after = Pt(6)

    headers = [
        "#",
        "Summary",
        "Length",
        "Penalty-Score",
        "Chrom / Gene & Transcript",
        "ForPos",
        "ForEnd",
        "RevPos",
        "RevEnd",
        "Seq (product)",
    ]
    ncols = len(headers)
    table = doc.add_table(rows=1 + len(amplicons), cols=ncols)
    table.style = "Table Grid"

    hdr_cells = table.rows[0].cells
    for j, title in enumerate(headers):
        _set_cell_run(hdr_cells[j], title, bold=True, size_pt=9)

    for i, amp in enumerate(amplicons):
        row_cells = table.rows[i + 1].cells
        values = [
            str(i + 1),
            extract_amplicon_summary(amp),
            str(amp.get("Length", "") or ""),
            format_penalty_score(amp.get("Penalty")),
            amplicon_chrom_label(amp),
            str(amp.get("ForPos", "") or ""),
            str(amp.get("ForEnd", "") or ""),
            str(amp.get("RevPos", "") or ""),
            str(amp.get("RevEnd", "") or ""),
            truncate_product_seq(amp.get("Seq")),
        ]
        for j, val in enumerate(values):
            _set_cell_run(row_cells[j], val, size_pt=8)

    doc.add_paragraph()


def create_primer_report(
    designResultsSummary_obj: DesignResultsSummary, selected_primer_index: int
):
    def create_var_annotation_block(var_info: AllelicVariantInfo) -> list[str]:
        annotation = []
        if var_info.gene_ID:
            annotation += [("Gene symbol", var_info.gene_symbol)]
            annotation += [("Gene-ID", var_info.gene_ID)]
        if isinstance(var_info, TranscriptVariantInfo):
            annotation += [("Transcript-ID", var_info.transcript_id)]
        annotation += [("HGVS", create_hgvs_notation(var_info))]
        return annotation

    prim_settings = designResultsSummary_obj.primer_settings
    var_info = designResultsSummary_obj.get_variant_info()
    primer_search_results = designResultsSummary_obj.get_primer_search_results()
    primer_pair = primer_search_results.primer_pairs[selected_primer_index - 1]

    # create new document
    doc = Document()
    # Layout definition
    section = doc.sections[0]
    section.top_margin = Inches(0.7)  # Oberer Rand
    section.bottom_margin = Inches(0.7)  # Unterer Rand
    section.left_margin = Inches(0.7)  # Linker Rand
    section.right_margin = Inches(0.7)  # Rechter Rand
    style = doc.styles["Normal"]
    font = style.font
    font.name = "Calibri"
    font.size = Pt(11)  # oder deine Wunschgröße, z. B. 12 pt

    # Add creation date to header, right-aligned
    header = section.header
    paragraph = header.paragraphs[0] if header.paragraphs else header.add_paragraph()
    paragraph.text = f"Created: {datetime.now().strftime('%d.%m.%Y %H:%M')}"
    paragraph.alignment = 2  # 2 = RIGHT

    doc.add_heading("Primer Designer Results", level=0)

    # add a paragraph with the variant information
    doc.add_heading("Input Variant:", level=1)
    for info, content in create_var_annotation_block(var_info):
        p = doc.add_paragraph(style="List Bullet")
        run = p.add_run(f"{info}: ")
        bold_run = p.add_run(content)
        bold_run.bold = True

    primer_data = (
        ("Sequence", primer_pair.left_seq, primer_pair.right_seq),
        (
            "(Relative) start, end",
            (f"{primer_pair.left_relPos_start}, {primer_pair.left_relPos_end}"),
            (f"{primer_pair.right_relPos_start}, {primer_pair.right_relPos_end}"),
        ),
        ("Tm", f"{primer_pair.tm[0]} °C", f"{primer_pair.tm[1]} °C"),
        ("GC-content", f"{primer_pair.gc[0]} %", f"{primer_pair.gc[1]} %"),
    )

    product_data = (
        ("Product size", f"{primer_pair.product_size} bp"),
        ("Product tm", f"{primer_pair.product_tm} °C "),
    )

    from primer_designer_app.utils.primer_utils import (
        INSILICO_ERROR,
        INSILICO_NOT_APPLICABLE,
        INSILICO_OK,
        INSILICO_OK_EMPTY,
    )

    def _insilico_report_line() -> str:
        st = getattr(primer_pair, "insilico_status", None)
        if st == INSILICO_NOT_APPLICABLE:
            return "Not applicable"
        if st == INSILICO_ERROR:
            return "Error (in-silico)"
        if st == INSILICO_OK_EMPTY:
            return "0"
        if st == INSILICO_OK:
            return str(len(primer_pair.amplicons or []))
        return str(len(primer_pair.amplicons or []))

    if getattr(prim_settings, "do_insilico_pcr", False):
        if prim_settings.context == "genomic":
            product_data = product_data + (
                ("Nr of amplicons (in Genome)", _insilico_report_line()),
            )
        elif prim_settings.context == "transcriptomic":
            product_data = product_data + (
                ("Nr of amplicons (in Transcriptome)", _insilico_report_line()),
            )

    # Creating a table object
    doc.add_heading("Primer Selection:", level=1)
    p = doc.add_paragraph("Link to all primer-pairs: ")
    webAppHost = settings.WEB_APP_HOST  # Load WEB_APP_HOST from settings
    # Load primer overview URL from urls.py
    primer_overview_url = reverse(
        "primer_designer_app:snv_indel_primers_overview_with_uuid",
        kwargs={"uuid": designResultsSummary_obj.id},
    )
    add_hyperlink(
        p,
        f"{webAppHost}{primer_overview_url}",
        "Primer results overview",
    )
    table = doc.add_table(rows=1, cols=3)

    # Adding heading in the 1st row of the table
    row = table.rows[0].cells
    row[1].text = ""
    paragraph = row[1].paragraphs[0]
    run = paragraph.add_run("Forward")
    run.bold = True

    paragraph = row[2].paragraphs[0]
    run = paragraph.add_run("Reverse")
    run.bold = True

    for criteria, Info_F, Info_R in primer_data:
        row = table.add_row().cells
        row[0].text = criteria
        row[1].text = Info_F
        row[2].text = Info_R

    for criteria, info in product_data:
        doc.add_paragraph(f"{criteria}: {info}", style="List Bullet")

    add_amplicon_detail_table_to_doc(doc, primer_pair, prim_settings)

    snp_analysis = getattr(designResultsSummary_obj, "snp_analysis_data", None) or {}
    if snp_analysis.get("enabled"):
        doc.add_heading("Known SNPs (gnomAD, MAF > 1%)", level=1)
        doc.add_paragraph(snp_analysis.get("message", ""))
        if snp_analysis.get("region"):
            region = snp_analysis["region"]
            doc.add_paragraph(
                f"Design region: chr{region.get('chromosome')}:"
                f"{region.get('start')}-{region.get('end')}",
                style="List Bullet",
            )
        pair_status = getattr(primer_pair, "snp_status", None)
        if pair_status:
            doc.add_paragraph(
                f"Selected primer pair SNP risk: {pair_status}",
                style="List Bullet",
            )
        conflicts = getattr(primer_pair, "snp_conflicts", None) or []
        if conflicts:
            doc.add_heading("SNP overlap with selected primers", level=2)
            table = doc.add_table(rows=1 + len(conflicts), cols=5)
            table.style = "Table Grid"
            for j, title in enumerate(["ID", "Genomic", "Alleles", "MAF", "Primer"]):
                _set_cell_run(table.rows[0].cells[j], title, bold=True, size_pt=9)
            for i, hit in enumerate(conflicts):
                row_cells = table.rows[i + 1].cells
                maf_val = hit.get("maf")
                maf_str = f"{maf_val:.4f}" if maf_val is not None else ""
                values = [
                    hit.get("id", ""),
                    f"{hit.get('genomic_start')}–{hit.get('genomic_end')}",
                    hit.get("alleles", ""),
                    maf_str,
                    hit.get("primer", ""),
                ]
                for j, val in enumerate(values):
                    _set_cell_run(row_cells[j], val, size_pt=8)

    # 3. Sequence visualization
    doc.add_heading("Sequence snippet:", level=1)
    doc.add_paragraph(
        "Sequence snippet is soft-masked: UTRs and introns are shown in lowercase "
        "letters, exons in uppercase letters."
    )

    doc.add_heading("Legend:", level=2)
    legend_content = [["Primer", 3], ["Target-Region", 5], ["Variant", 4]]

    for region, color in legend_content:
        list_obj = doc.add_paragraph(style="List Bullet")
        run = list_obj.add_run(region)
        run.font.highlight_color = color

    paragraph = doc.add_paragraph()
    visualize_sequence_as_docx(paragraph, prim_settings, var_info, primer_pair)

    # Save the document to an in-memory file
    buffer = io.BytesIO()
    doc.save(buffer)
    buffer.seek(0)
    return buffer


def _init_report_document(title: str) -> Document:
    doc = Document()
    section = doc.sections[0]
    section.top_margin = Inches(0.7)
    section.bottom_margin = Inches(0.7)
    section.left_margin = Inches(0.7)
    section.right_margin = Inches(0.7)
    style = doc.styles["Normal"]
    style.font.name = "Calibri"
    style.font.size = Pt(11)

    header = section.header
    paragraph = header.paragraphs[0] if header.paragraphs else header.add_paragraph()
    paragraph.text = f"Created: {datetime.now().strftime('%d.%m.%Y %H:%M')}"
    paragraph.alignment = 2

    doc.add_heading(title, level=0)
    return doc


def _add_bullet_field(doc: Document, label: str, value: str) -> None:
    p = doc.add_paragraph(style="List Bullet")
    p.add_run(f"{label}: ")
    bold_run = p.add_run(value)
    bold_run.bold = True


def _add_sv_primer_pairs_table(doc: Document, primer_rows: list) -> None:
    headers = [
        "Rank",
        "Forward primer",
        "Reverse primer",
        "Penalty",
        "Product size (bp)",
        "Primer Tm (°C)",
        "Primer GC (%)",
        "Forward genomic",
        "Reverse genomic",
    ]
    if not primer_rows:
        doc.add_paragraph("No primer pairs found for this window.")
        return

    table = doc.add_table(rows=1 + len(primer_rows), cols=len(headers))
    table.style = "Table Grid"
    for j, title in enumerate(headers):
        _set_cell_run(table.rows[0].cells[j], title, bold=True, size_pt=9)

    for i, row in enumerate(primer_rows):
        pair = row["pair"]
        genomic = row["genomic_positions"]
        tm_text = (
            f"{pair.tm[0]}, {pair.tm[1]}" if pair.tm and len(pair.tm) >= 2 else "—"
        )
        gc_text = (
            f"{pair.gc[0]}, {pair.gc[1]}" if pair.gc and len(pair.gc) >= 2 else "—"
        )
        values = [
            str(i + 1),
            pair.left_seq,
            pair.right_seq,
            str(pair.penalty),
            str(pair.product_size),
            tm_text,
            gc_text,
            f"{genomic['forward_start']}–{genomic['forward_end']}",
            f"{genomic['reverse_start']}–{genomic['reverse_end']}",
        ]
        for j, val in enumerate(values):
            _set_cell_run(table.rows[i + 1].cells[j], val, size_pt=8)


def create_structural_variant_primer_report(
    design_results_summary: DesignResultsSummary,
) -> io.BytesIO:
    """
    Word report for structural variant primer design: all primer pairs per window,
    without sequence visualization or primer-pair selection.
    """
    sv_info = design_results_summary.get_structural_variant_info_data()
    sv_results = design_results_summary.get_sv_primer_results()
    prim_settings = design_results_summary.primer_settings

    doc = _init_report_document("Structural Variant Primer Designer Results")

    doc.add_heading("Structural variant query", level=1)
    chromosome = sv_info.get("chromosome", "")
    start_pos = sv_info.get("start_position", "")
    end_pos = sv_info.get("end_position", "")
    _add_bullet_field(doc, "Reference genome", sv_info.get("reference_genome", ""))
    _add_bullet_field(doc, "Chromosome", f"chr{chromosome}")
    _add_bullet_field(doc, "Start position", str(start_pos))
    _add_bullet_field(doc, "End position", str(end_pos))
    if start_pos and end_pos:
        span = int(end_pos) - int(start_pos) + 1
        _add_bullet_field(doc, "Span", f"{span} bp")

    doc.add_heading("Design windows", level=2)
    for window in sv_info.get("windows", []):
        doc.add_paragraph(
            f"{window.get('label', '').replace('_', ' ').title()}: "
            f"{window.get('window_start_genomic')}–{window.get('window_end_genomic')} "
            f"(chr{chromosome})",
            style="List Bullet",
        )

    doc.add_heading("Primer design parameters", level=1)
    if prim_settings:
        _add_bullet_field(doc, "Optimal Tm", f"{prim_settings.tm} °C")
        _add_bullet_field(doc, "GC target", f"{prim_settings.gc} %")
        _add_bullet_field(doc, "Max poly-X", str(prim_settings.max_poly_x))
        if prim_settings.productsize_range:
            _add_bullet_field(
                doc,
                "Product size range",
                f"{prim_settings.productsize_range[0]}–{prim_settings.productsize_range[1]} bp",
            )

    doc.add_heading("Designed primer pairs", level=1)
    doc.add_paragraph(
        "All primer pairs returned by Primer3 for each design window are listed below. "
        "Genomic coordinates refer to the selected reference assembly."
    )

    labels_in_report = [label for label in SV_WINDOW_ORDER if label in sv_results]
    labels_in_report.extend(
        label for label in sv_results if label not in labels_in_report
    )

    for label in labels_in_report:
        window_result = sv_results[label]
        window = window_result["window"]
        window_title = label.replace("_", " ").title()
        doc.add_heading(window_title, level=2)
        doc.add_paragraph(
            f"Window coordinates (chr{chromosome}): "
            f"{window['window_start_genomic']}–{window['window_end_genomic']}"
        )
        _add_sv_primer_pairs_table(doc, window_result.get("primer_rows", []))
        doc.add_paragraph()

    buffer = io.BytesIO()
    doc.save(buffer)
    buffer.seek(0)
    return buffer
