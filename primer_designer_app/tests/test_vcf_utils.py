import io
import unittest

from primer_designer_app.utils.vcf_utils import (
    VcfRecord,
    compute_fetch_window,
    parse_vcf_upload,
    spike_vcf_variants,
    template_range_for_genomic,
)


class VcfUtilsTests(unittest.TestCase):
    def test_parse_vcf_filters_chromosome(self):
        vcf = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT
7\t100\t.\tA\tG
7\t200\trs2\tC\tT
8\t50\t.\tG\tA
"""
        records = parse_vcf_upload(io.BytesIO(vcf.encode()), "7")
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].pos, 100)

    def test_spike_snv_and_indel(self):
        ref = "AAAA" + "CCCC" + "TTTT"  # positions 1-12
        records = [
            VcfRecord(chrom="7", pos=5, ref="C", alt="G", rsid="rsA"),
            VcfRecord(chrom="7", pos=10, ref="T", alt="TA", rsid="rsB"),
        ]
        spiked, applied, deltas = spike_vcf_variants(ref, 1, records)
        self.assertEqual(len(applied), 2)
        self.assertEqual(spiked[4], "G")
        self.assertIn("TA", spiked)

    def test_skip_primary_interval(self):
        ref = "ACGTACGT"
        records = [VcfRecord(chrom="7", pos=3, ref="G", alt="A", rsid="rs1")]
        spiked, applied, _ = spike_vcf_variants(ref, 1, records, skip_interval=(3, 3))
        self.assertEqual(spiked, ref)
        self.assertEqual(applied, [])

    def test_template_range_after_spike(self):
        region_start = 100
        primary_pos = 105
        records = [
            VcfRecord(chrom="7", pos=101, ref="A", alt="AA", rsid="ins"),
        ]
        ref = "A" * 20
        _, _, deltas = spike_vcf_variants(ref, region_start, records)
        t_start, t_end = template_range_for_genomic(
            region_start, primary_pos, primary_pos, deltas
        )
        self.assertEqual(t_start, t_end)
        self.assertEqual(t_start, 5 + 1)  # pos 105 -> index 5, +1 from ins before

    def test_template_range_snv_not_inverted_after_upstream_indel(self):
        """SNV primary interval must not collapse to empty ref_bases (false INS)."""
        region_start = 100
        primary_pos = 105
        records = [VcfRecord(chrom="7", pos=101, ref="A", alt="AA", rsid="ins")]
        ref = "A" * 20
        _, _, deltas = spike_vcf_variants(ref, region_start, records)
        t_start, t_end = template_range_for_genomic(
            region_start, primary_pos, primary_pos, deltas
        )
        self.assertGreaterEqual(t_end, t_start)

    def test_compute_fetch_window(self):
        start, end = compute_fetch_window(
            1000,
            1000,
            [VcfRecord("7", 800, "A", "G"), VcfRecord("7", 1200, "C", "T")],
            flank=100,
        )
        self.assertEqual(start, 700)
        self.assertEqual(end, 1300)


if __name__ == "__main__":
    unittest.main()
