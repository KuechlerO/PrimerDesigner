import unittest

from types import SimpleNamespace

from primer_designer_app.utils.display_utils import (
    DISPLAY_FLANK,
    REPORT_DISPLAY_FLANK,
    compute_display_bounds,
    compute_report_display_bounds,
)


def _pair(f0, f1, r0, r1):
    return SimpleNamespace(
        left_relPos_start=f0,
        left_relPos_end=f1,
        right_relPos_start=r0,
        right_relPos_end=r1,
    )


class DisplayBoundsTests(unittest.TestCase):
    def test_window_spans_both_primers(self):
        pair = _pair(80, 99, 1700, 1719)
        start, end = compute_display_bounds(2001, 1000, 1000, 950, 101, [pair])
        self.assertLessEqual(start, 80)
        self.assertGreaterEqual(end, 1720)
        self.assertLess(end - start, 2001)

    def test_js_coords_after_offset(self):
        pair = _pair(80050, 80069, 80150, 80169)
        offset, _ = compute_display_bounds(120000, 80080, 80080, 80030, 101, [pair])
        self.assertGreaterEqual(80050 - offset, 0)
        self.assertLess(80169 - offset, 120000 - offset)

    def test_report_flank_wider_than_ui(self):
        pair = _pair(800, 819, 1500, 1519)
        ui_start, ui_end = compute_display_bounds(5000, 1000, 1000, 950, 101, [pair])
        rep_start, rep_end = compute_report_display_bounds(
            5000, 1000, 1000, 950, 101, pair
        )
        self.assertLess(rep_start, ui_start)
        self.assertGreater(rep_end, ui_end)
        self.assertGreaterEqual(
            ui_start - rep_start, REPORT_DISPLAY_FLANK - DISPLAY_FLANK
        )


if __name__ == "__main__":
    unittest.main()
