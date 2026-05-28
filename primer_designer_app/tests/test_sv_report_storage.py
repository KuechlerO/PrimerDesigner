import types
import unittest

from primer_designer_app.utils.sv_storage import (
    SV_DESIGN_TYPE,
    deserialize_sv_results_from_storage,
    serialize_structural_variant_info,
    serialize_sv_results_for_storage,
)


class SvReportStorageTests(unittest.TestCase):
    def test_roundtrip_serialization(self):
        class _FakeInfo:
            chromosome = "12"
            start_position = 100
            end_position = 200
            reference_genome = "GRCh37"
            windows = [
                types.SimpleNamespace(
                    label="upstream",
                    window_start_genomic=1,
                    window_end_genomic=99,
                )
            ]

        class _FakePair:
            def to_dict(self):
                return {
                    "index": 0,
                    "left_seq": "AAA",
                    "right_seq": "TTT",
                    "penalty": 0.2,
                    "product_size": 120,
                    "product_tm": 60.0,
                    "left_relPos_start": 1,
                    "left_relPos_end": 20,
                    "right_relPos_start": 80,
                    "right_relPos_end": 99,
                    "gc": [50.0, 51.0],
                    "tm": [60.0, 59.0],
                    "insilico_seq": "",
                    "amplicons": [],
                    "insilico_status": "not_applicable",
                    "insilico_error_detail": None,
                }

        results = {
            "upstream": {
                "design_window": _FakeInfo.windows[0],
                "primer_rows": [
                    {
                        "pair": _FakePair(),
                        "genomic_positions": {
                            "forward_start": 10,
                            "forward_end": 29,
                            "reverse_start": 80,
                            "reverse_end": 99,
                        },
                    }
                ],
            }
        }

        stored_info = serialize_structural_variant_info(_FakeInfo)
        self.assertEqual(stored_info["design_type"], SV_DESIGN_TYPE)

        stored_results = serialize_sv_results_for_storage(results)
        payload = {"design_type": SV_DESIGN_TYPE, "windows": stored_results}
        restored = deserialize_sv_results_from_storage(payload, lambda d: d["left_seq"])
        self.assertIn("upstream", restored)
        self.assertEqual(restored["upstream"]["primer_rows"][0]["pair"], "AAA")


if __name__ == "__main__":
    unittest.main()
