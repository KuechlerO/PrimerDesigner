import unittest

from types import SimpleNamespace

from primer_designer_app.utils.display_utils import (
    DISPLAY_FLANK,
    compute_display_bounds,
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


if __name__ == "__main__":
    unittest.main()
