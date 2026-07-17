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

    def test_parse_vcf_skips_structural_alleles(self):
        vcf = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT
7\t100\t.\tA\t<DEL>
7\t200\t.\t<INS>\tA
7\t300\trs5\tG\tT
7\t400\t.\tG\t*
"""
        records = parse_vcf_upload(io.BytesIO(vcf.encode()), "7")
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].pos, 300)
        self.assertEqual(records[0].rsid, "rs5")

    def test_parse_vcf_skips_non_nucleotide_alleles(self):
        vcf = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT
7\t100\t.\tACGTRYK\tA
7\t200\t.\tA\tACGTRYK
7\t300\t.\tA\tG
"""
        records = parse_vcf_upload(io.BytesIO(vcf.encode()), "7")
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].pos, 300)

    def test_spike_skips_on_ref_mismatch(self):
        ref = "AAAACCCCTTTT"
        records = [VcfRecord(chrom="7", pos=5, ref="G", alt="T", rsid="rsMismatch")]
        spiked, applied, deltas = spike_vcf_variants(ref, 1, records)
        self.assertEqual(spiked, ref)
        self.assertEqual(applied, [])
        self.assertEqual(deltas, [])

    def test_spike_skips_outside_fetched_window(self):
        ref = "AAAA"
        records = [VcfRecord(chrom="7", pos=100, ref="A", alt="G", rsid="rsFar")]
        spiked, applied, _ = spike_vcf_variants(ref, 1, records)
        self.assertEqual(spiked, ref)
        self.assertEqual(applied, [])

    def test_multiple_indels_compound_offset(self):
        # Two upstream insertions before a downstream SNV: template offset for the
        # SNV must account for both insertions' net length changes.
        ref = "A" * 5 + "C" * 5 + "T" * 5  # positions 1-15
        records = [
            VcfRecord(chrom="7", pos=2, ref="A", alt="AAA", rsid="ins1"),  # +2
            VcfRecord(chrom="7", pos=6, ref="C", alt="CC", rsid="ins2"),  # +1
            VcfRecord(chrom="7", pos=11, ref="T", alt="G", rsid="snv"),
        ]
        spiked, applied, deltas = spike_vcf_variants(ref, 1, records)
        self.assertEqual(len(applied), 3)
        t_start, t_end = template_range_for_genomic(1, 11, 11, deltas)
        # Base index for pos 11 (0-based, no spikes) would be 10; +3 from upstream
        # insertions (ins1 +2, ins2 +1) => 13.
        self.assertEqual(t_start, 13)
        self.assertEqual(t_end, 13)
        self.assertEqual(spiked[13], "G")

    def test_normalize_chromosome_variants(self):
        from primer_designer_app.utils.vcf_utils import normalize_chromosome

        self.assertEqual(normalize_chromosome("chr7"), "7")
        self.assertEqual(normalize_chromosome("CHRX"), "X")
        self.assertEqual(normalize_chromosome("MT"), "M")
        self.assertEqual(normalize_chromosome(" chr1 "), "1")
        self.assertEqual(normalize_chromosome(""), "")


if __name__ == "__main__":
    unittest.main()
