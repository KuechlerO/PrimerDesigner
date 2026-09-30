import unittest
from unittest.mock import MagicMock, patch

from primer_designer_app.utils.variant_info import (
    AllelicVariantInfo,
    IndelType,
    SequenceVariantInfo,
    StructuralVariantInfo,
    StructuralVariantWindow,
)


class StructuralVariantInfoCreateDesignWindowsTests(unittest.TestCase):
    def test_typical_deletion(self):
        sv_info = StructuralVariantInfo(
            chromosome="1",
            start_position=1000,
            end_position=5000,
            reference_genome="GRCh37",
        )
        windows = sv_info.create_design_windows()
        by_label = {w.label: w for w in windows}

        self.assertEqual(
            set(by_label), {"upstream", "downstream", "internal_1", "internal_2"}
        )
        # upstream: [start - flank, start - 1], flank default 5000
        self.assertEqual(by_label["upstream"].window_start_genomic, 1)
        self.assertEqual(by_label["upstream"].window_end_genomic, 999)
        # downstream: [end + 1, end + flank]
        self.assertEqual(by_label["downstream"].window_start_genomic, 5001)
        self.assertEqual(by_label["downstream"].window_end_genomic, 10000)
        # internal windows: min(2000, sv_length/2) = min(2000, 2000.5) = 2000
        self.assertEqual(by_label["internal_1"].window_start_genomic, 1000)
        self.assertEqual(by_label["internal_1"].window_end_genomic, 2999)
        self.assertEqual(by_label["internal_2"].window_start_genomic, 3001)
        self.assertEqual(by_label["internal_2"].window_end_genomic, 5000)

    def test_sv_length_50_boundary_succeeds(self):
        # length exactly 50 -> half_sv_length == 25 -> effective_internal_window_size == 25
        # (not < 25), so this must succeed.
        sv_info = StructuralVariantInfo(
            chromosome="1",
            start_position=100,
            end_position=149,
            reference_genome="GRCh37",
        )
        self.assertEqual(sv_info.structural_variant_length, 50)
        windows = sv_info.create_design_windows()
        self.assertEqual(len(windows), 4)

    def test_sv_length_below_50_raises(self):
        sv_info = StructuralVariantInfo(
            chromosome="1",
            start_position=100,
            end_position=148,  # length 49
            reference_genome="GRCh37",
        )
        with self.assertRaisesRegex(
            ValueError, "Structural variant must span at least 50 bases"
        ):
            sv_info.create_design_windows()

    def test_upstream_flank_clamped_to_chromosome_start(self):
        # start_position - flank_window_size goes negative -> clamped to 1.
        sv_info = StructuralVariantInfo(
            chromosome="1",
            start_position=2000,
            end_position=6000,
            reference_genome="GRCh37",
        )
        windows = sv_info.create_design_windows(flank_window_size=5000)
        upstream = next(w for w in windows if w.label == "upstream")
        self.assertEqual(upstream.window_start_genomic, 1)
        self.assertEqual(upstream.window_end_genomic, 1999)

    def test_start_position_at_chromosome_start_raises(self):
        # start_position == 1 -> upstream_window_end == 0 -> cannot create upstream window.
        sv_info = StructuralVariantInfo(
            chromosome="1",
            start_position=1,
            end_position=100,
            reference_genome="GRCh37",
        )
        with self.assertRaisesRegex(
            ValueError, "Upstream window cannot be created at chromosome start"
        ):
            sv_info.create_design_windows()

    def test_internal_window_too_small_raises(self):
        sv_info = StructuralVariantInfo(
            chromosome="1",
            start_position=1000,
            end_position=6000,  # length 5001, half=2500.5
            reference_genome="GRCh37",
        )
        with self.assertRaisesRegex(
            ValueError, "Internal window size is too small for primer design"
        ):
            sv_info.create_design_windows(internal_window_size=10)

    def test_start_greater_than_end_raises(self):
        with self.assertRaisesRegex(
            ValueError, "Start position must not exceed end position"
        ):
            StructuralVariantInfo(
                chromosome="1",
                start_position=200,
                end_position=100,
                reference_genome="GRCh37",
            )

    def test_non_positive_positions_raise(self):
        with self.assertRaisesRegex(
            ValueError, "Genomic positions must be positive integers"
        ):
            StructuralVariantInfo(
                chromosome="1",
                start_position=0,
                end_position=100,
                reference_genome="GRCh37",
            )


class StructuralVariantWindowTargetTests(unittest.TestCase):
    def _window(self) -> StructuralVariantWindow:
        return StructuralVariantWindow(
            label="upstream",
            window_start_genomic=1,
            window_end_genomic=1000,  # window_length == 1000
        )

    def test_set_target_valid(self):
        window = self._window()
        window.set_target(100, 200)
        self.assertEqual(window.target_start_in_window, 100)
        self.assertEqual(window.target_length, 200)
        self.assertEqual(window.get_primer3_target(), [100, 200])

    def test_set_target_negative_start_raises(self):
        window = self._window()
        with self.assertRaisesRegex(ValueError, "target_start_in_window must be >= 0"):
            window.set_target(-1, 100)

    def test_set_target_non_positive_length_raises(self):
        window = self._window()
        with self.assertRaisesRegex(ValueError, "target_length must be > 0"):
            window.set_target(0, 0)

    def test_set_target_exceeds_window_bounds_raises(self):
        window = self._window()
        with self.assertRaisesRegex(ValueError, "Target exceeds window boundaries"):
            window.set_target(900, 200)

    def test_set_default_target_centers_in_window(self):
        window = self._window()
        window.set_default_target(target_length=200)
        # (1000 - 200) / 2 == 400
        self.assertEqual(window.target_start_in_window, 400)
        self.assertEqual(window.target_length, 200)

    def test_set_default_target_clamps_to_window_length(self):
        window = StructuralVariantWindow(
            label="internal_1",
            window_start_genomic=1,
            window_end_genomic=50,  # window_length == 50
        )
        window.set_default_target(target_length=150)
        self.assertEqual(window.target_length, 50)
        self.assertEqual(window.target_start_in_window, 0)

    def test_get_primer3_target_before_set_raises(self):
        window = self._window()
        with self.assertRaisesRegex(ValueError, "Primer3 target has not been set"):
            window.get_primer3_target()


class SequenceVariantInfoParseTests(unittest.TestCase):
    def test_snv_annotation(self):
        var = SequenceVariantInfo(input_seq="AC[G>A]GT", ref_genome="GRCh37")
        self.assertEqual(var.ref_seq, "ACGGT")
        self.assertEqual(var.ref_bases, "G")
        self.assertEqual(var.new_bases, "A")
        self.assertEqual(var.relative_pos, (2, 2))
        self.assertEqual(var.indel_type, IndelType.SNV)

    def test_insertion_annotation(self):
        var = SequenceVariantInfo(input_seq="AC[-/GG]GT", ref_genome="GRCh37")
        self.assertEqual(var.ref_seq, "ACGT")
        self.assertEqual(var.ref_bases, "")
        self.assertEqual(var.new_bases, "GG")
        self.assertEqual(var.relative_pos, (2, 2))
        self.assertEqual(var.indel_type, IndelType.INS)

    def test_deletion_annotation(self):
        var = SequenceVariantInfo(input_seq="AC[GA/-]GT", ref_genome="GRCh37")
        self.assertEqual(var.ref_seq, "ACGAGT")
        self.assertEqual(var.ref_bases, "GA")
        self.assertEqual(var.new_bases, "")
        self.assertEqual(var.relative_pos, (2, 3))
        self.assertEqual(var.indel_type, IndelType.DEL)

    def test_delins_annotation(self):
        var = SequenceVariantInfo(input_seq="AC[GA/CT]GT", ref_genome="GRCh37")
        self.assertEqual(var.ref_seq, "ACGAGT")
        self.assertEqual(var.ref_bases, "GA")
        self.assertEqual(var.new_bases, "CT")
        self.assertEqual(var.relative_pos, (2, 3))
        self.assertEqual(var.indel_type, IndelType.DELINS)

    def test_missing_brackets_raises(self):
        with self.assertRaisesRegex(ValueError, "Sequence annotation missing brackets"):
            SequenceVariantInfo(input_seq="ACGT", ref_genome="GRCh37")

    def test_strips_whitespace_and_newlines(self):
        var = SequenceVariantInfo(input_seq="AC\n[G>A]\nGT", ref_genome="GRCh37")
        self.assertEqual(var.ref_seq, "ACGGT")
        self.assertEqual(var.indel_type, IndelType.SNV)


class AllelicVariantInfoTests(unittest.TestCase):
    def test_requires_relative_pos(self):
        with self.assertRaisesRegex(
            ValueError, "relative_pos must be provided for AllelicVariantInfo"
        ):
            AllelicVariantInfo(ref_seq="AAAA")

    def test_determine_indel_type_snv(self):
        var = AllelicVariantInfo(
            ref_seq="AAAGAAA", ref_bases="G", new_bases="A", relative_pos=(3, 3)
        )
        self.assertEqual(var._determine_indel_type(), IndelType.SNV)

    def test_determine_indel_type_ins(self):
        var = AllelicVariantInfo(
            ref_seq="AAAAAA", ref_bases="", new_bases="GG", relative_pos=(3, 3)
        )
        self.assertEqual(var._determine_indel_type(), IndelType.INS)

    def test_determine_indel_type_del(self):
        var = AllelicVariantInfo(
            ref_seq="AAAGAAAA", ref_bases="GA", new_bases="", relative_pos=(3, 4)
        )
        self.assertEqual(var._determine_indel_type(), IndelType.DEL)

    def test_determine_indel_type_delins(self):
        var = AllelicVariantInfo(
            ref_seq="AAAGAAAA", ref_bases="GA", new_bases="CT", relative_pos=(3, 4)
        )
        self.assertEqual(var._determine_indel_type(), IndelType.DELINS)

    def test_get_seq_mutated_and_input_snv(self):
        var = AllelicVariantInfo(
            ref_seq="AAAGAAA",
            ref_bases="G",
            new_bases="A",
            relative_pos=(3, 3),
            indel_type=IndelType.SNV,
        )
        self.assertEqual(var.get_seq("mutated"), "AAAAAAA")
        self.assertEqual(var.get_seq("input"), "AAA[G>A]AAA")

    def test_get_seq_mutated_and_input_ins(self):
        var = AllelicVariantInfo(
            ref_seq="AAAAAA",
            ref_bases="",
            new_bases="GG",
            relative_pos=(3, 3),
            indel_type=IndelType.INS,
        )
        self.assertEqual(var.get_seq("mutated"), "AAAGGAAA")
        self.assertEqual(var.get_seq("input"), "AAA[-/GG]AAA")

    def test_get_seq_mutated_and_input_del(self):
        var = AllelicVariantInfo(
            ref_seq="AAAGAAAA",
            ref_bases="GA",
            new_bases="",
            relative_pos=(3, 4),
            indel_type=IndelType.DEL,
        )
        self.assertEqual(var.get_seq("mutated"), "AAAAAA")
        self.assertEqual(var.get_seq("input"), "AAA[GA/-]AAA")

    def test_get_seq_mutated_and_input_delins(self):
        var = AllelicVariantInfo(
            ref_seq="AAAGAAAA",
            ref_bases="GA",
            new_bases="CT",
            relative_pos=(3, 4),
            indel_type=IndelType.DELINS,
        )
        self.assertEqual(var.get_seq("mutated"), "AAACTAAA")
        self.assertEqual(var.get_seq("input"), "AAA[GA/CT]AAA")

    def test_get_seq_invalid_output_type_raises(self):
        var = AllelicVariantInfo(
            ref_seq="AAAGAAA",
            ref_bases="G",
            new_bases="A",
            relative_pos=(3, 3),
            indel_type=IndelType.SNV,
        )
        with self.assertRaisesRegex(ValueError, "Invalid output type"):
            var.get_seq("bogus")

    def test_get_genomic_pos(self):
        var = AllelicVariantInfo(
            ref_seq="AAAA",
            relative_pos=(0, 0),
            genomic_pos={"chr": "7", "pos": [100, 100]},
        )
        self.assertEqual(var.get_genomic_pos(), [100, 100])

    def test_set_attribute(self):
        var = AllelicVariantInfo(ref_seq="AAAA", relative_pos=(0, 0))
        var.set_attribute("gene_symbol", "TP53")
        self.assertEqual(var.gene_symbol, "TP53")

    def test_ignores_unknown_kwargs(self):
        # Should not raise even though AS_seq is not a dataclass field.
        var = AllelicVariantInfo(ref_seq="AAAA", relative_pos=(0, 0), AS_seq="ignored")
        self.assertFalse(hasattr(var, "AS_seq"))


class CdsOffsetOnCdnaTests(unittest.TestCase):
    def test_finds_cds_within_cdna(self):
        from primer_designer_app.utils.variant_info import cds_offset_on_cdna

        cdna = "aaaatttATGCGTAAAggg"
        cds = "ATGCGTAAA"
        self.assertEqual(cds_offset_on_cdna(cdna, cds), 7)

    def test_case_insensitive(self):
        from primer_designer_app.utils.variant_info import cds_offset_on_cdna

        self.assertEqual(cds_offset_on_cdna("nnnATG", "atg"), 3)

    def test_missing_cds_raises(self):
        from primer_designer_app.utils.variant_info import cds_offset_on_cdna

        with self.assertRaisesRegex(ValueError, "Could not locate CDS"):
            cds_offset_on_cdna("AAAA", "TTTT")


class TranscriptVariantInfoCdsRemapTests(unittest.TestCase):
    def _mock_client(self):
        client = MagicMock()
        client.get_gene_symbol_for_transcriptID.return_value = (
            "PPAPDC1B",
            "ENSG00000147535",
            "ENST00000424479.2",
        )
        # 5' UTR length 100, then CDS
        cds = "ATG" + ("C" * 200) + "TAA"
        cdna = ("a" * 100) + cds + ("g" * 50)
        client.get_transcript_sequence.side_effect = (
            lambda tid, seq_type, mask_feature=None: (
                cdna if seq_type == "cdna" else cds
            )
        )
        client.map_coordinates.return_value = {
            "mappings": [
                {
                    "start": 1000,
                    "end": 1000,
                    "seq_region_name": "8",
                    "strand": 1,
                }
            ]
        }
        client.get_overlapped_genes_details_for_region.return_value = {
            "ENSG1": "PPAPDC1B"
        }
        return client

    @patch("primer_designer_app.utils.variant_info.EnsemblClient")
    def test_cds_input_remaps_onto_cdna_template(self, mock_cls):
        from primer_designer_app.utils.variant_info import (
            ReferenceType,
            TranscriptVariantInfo,
        )

        mock_cls.return_value = self._mock_client()
        # c.5 → 0-based CDS index 4
        var = TranscriptVariantInfo(
            transcript_id="ENST00000424479.2",
            ref_genome="GRCh37",
            reference_type=ReferenceType.CDS,
            new_bases="T",
            relative_pos=(4, 4),
        )
        self.assertEqual(var.input_relative_pos, (4, 4))
        self.assertEqual(var.relative_pos, (104, 104))
        self.assertEqual(len(var.ref_seq), 100 + 206 + 50)
        self.assertTrue(var.ref_seq.startswith("a" * 100))
        # Genomic map uses remapped cDNA coords (Ensembl /map/cdna/ is faster than /map/cds/)
        mock_cls.return_value.map_coordinates.assert_called_with(
            "ENST00000424479.2", 104, 104, "cdna"
        )

    @patch("primer_designer_app.utils.variant_info.EnsemblClient")
    def test_cds_remap_is_idempotent_on_reload(self, mock_cls):
        from primer_designer_app.utils.variant_info import (
            ReferenceType,
            TranscriptVariantInfo,
        )

        mock_cls.return_value = self._mock_client()
        first = TranscriptVariantInfo(
            transcript_id="ENST00000424479.2",
            ref_genome="GRCh37",
            reference_type=ReferenceType.CDS,
            new_bases="T",
            relative_pos=(4, 4),
        )
        # Simulate deserialize with stored fields (including remapped relative_pos).
        # Must NOT call Ensembl again when stored ref_seq/genomic_pos/ref_bases exist.
        mock_cls.return_value = self._mock_client()
        second = TranscriptVariantInfo(
            transcript_id=first.transcript_id,
            ref_genome=first.ref_genome,
            reference_type=ReferenceType.CDS,
            new_bases=first.new_bases,
            relative_pos=first.relative_pos,
            input_relative_pos=first.input_relative_pos,
            ref_bases=first.ref_bases,
            ref_seq=first.ref_seq,
            gene_ID=first.gene_ID,
            gene_symbol=first.gene_symbol,
            genomic_pos=first.genomic_pos,
            indel_type=first.indel_type,
        )
        self.assertEqual(second.input_relative_pos, (4, 4))
        self.assertEqual(second.relative_pos, (104, 104))
        self.assertEqual(second.relative_pos, first.relative_pos)
        mock_cls.return_value.get_transcript_sequence.assert_not_called()
        mock_cls.return_value.get_gene_symbol_for_transcriptID.assert_not_called()

    @patch("primer_designer_app.utils.variant_info.EnsemblClient")
    def test_cdna_mode_keeps_input_coords(self, mock_cls):
        from primer_designer_app.utils.variant_info import (
            ReferenceType,
            TranscriptVariantInfo,
        )

        mock_cls.return_value = self._mock_client()
        var = TranscriptVariantInfo(
            transcript_id="ENST00000424479.2",
            ref_genome="GRCh37",
            reference_type=ReferenceType.CDNA,
            new_bases="T",
            relative_pos=(50, 50),
        )
        self.assertEqual(var.input_relative_pos, (50, 50))
        self.assertEqual(var.relative_pos, (50, 50))
        mock_cls.return_value.map_coordinates.assert_called_with(
            "ENST00000424479.2", 50, 50, "cdna"
        )


class TranscriptHgvsTests(unittest.TestCase):
    @patch("primer_designer_app.utils.variant_info.EnsemblClient")
    def test_hgvs_uses_one_based_input_relative_pos(self, mock_cls):
        from primer_designer_app.utils.helpers import create_hgvs_notation
        from primer_designer_app.utils.variant_info import (
            ReferenceType,
            TranscriptVariantInfo,
        )

        client = MagicMock()
        client.get_gene_symbol_for_transcriptID.return_value = (
            "PPAPDC1B",
            "ENSG00000147535",
            "ENST00000424479.2",
        )
        cds = "ATGCG" + ("C" * 100)
        cdna = ("a" * 100) + cds
        client.get_transcript_sequence.side_effect = (
            lambda tid, seq_type, mask_feature=None: (
                cdna if seq_type == "cdna" else cds
            )
        )
        client.map_coordinates.return_value = {
            "mappings": [
                {
                    "start": 1000,
                    "end": 1000,
                    "seq_region_name": "8",
                    "strand": 1,
                }
            ]
        }
        client.get_overlapped_genes_details_for_region.return_value = {
            "ENSG1": "PPAPDC1B"
        }
        mock_cls.return_value = client

        var = TranscriptVariantInfo(
            transcript_id="ENST00000424479.2",
            ref_genome="GRCh37",
            reference_type=ReferenceType.CDS,
            new_bases="T",
            relative_pos=(4, 4),
        )
        hgvs = create_hgvs_notation(var)
        self.assertIn("c.5G>T", hgvs)
        self.assertIn("PPAPDC1B(ENST00000424479.2)", hgvs)


if __name__ == "__main__":
    unittest.main()
