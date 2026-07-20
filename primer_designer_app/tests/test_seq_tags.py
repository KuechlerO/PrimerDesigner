import unittest

from primer_designer_app.templatetags.seq_tags import chunk_html
from primer_designer_app.utils.display_utils import DEFAULT_CHUNK_WIDTH


class ChunkHtmlTests(unittest.TestCase):
    def test_plain_text_split_into_chunks(self):
        seq = "ACGT" * 10  # 40 chars
        chunks = chunk_html(seq, width=10)
        self.assertEqual(len(chunks), 4)
        self.assertEqual(chunks[0]["start"], 1)
        self.assertEqual(chunks[0]["chunk"], "ACGTACGTAC")
        self.assertEqual(chunks[1]["start"], 11)
        self.assertEqual(chunks[3]["start"], 31)

    def test_plain_text_shorter_than_width_single_chunk(self):
        chunks = chunk_html("ACGT", width=10)
        self.assertEqual(len(chunks), 1)
        self.assertEqual(chunks[0], {"start": 1, "chunk": "ACGT"})

    def test_tag_spanning_chunk_boundary_is_reopened(self):
        # 8 bases inside a <mark> tag, width=5 -> boundary falls inside the tag.
        seq_html = "AAA<mark>AAAAA</mark>AA"
        chunks = chunk_html(seq_html, width=5)
        self.assertEqual(len(chunks), 2)
        # First chunk should close the still-open <mark> tag.
        self.assertIn("<mark>", chunks[0]["chunk"])
        self.assertIn("</mark>", chunks[0]["chunk"])
        # Second chunk should reopen <mark> for the remaining highlighted bases.
        self.assertTrue(chunks[1]["chunk"].startswith("<mark>"))
        self.assertIn("</mark>", chunks[1]["chunk"])
        self.assertEqual(chunks[1]["start"], 6)

    def test_none_input_returns_empty_list(self):
        self.assertEqual(chunk_html(None), [])

    def test_width_zero_falls_back_to_default(self):
        seq = "A" * (DEFAULT_CHUNK_WIDTH + 5)
        chunks = chunk_html(seq, width=0)
        self.assertEqual(len(chunks), 2)
        self.assertEqual(len(chunks[0]["chunk"]), DEFAULT_CHUNK_WIDTH)

    def test_negative_width_falls_back_to_default(self):
        seq = "A" * (DEFAULT_CHUNK_WIDTH + 5)
        chunks = chunk_html(seq, width=-1)
        self.assertEqual(len(chunks), 2)
        self.assertEqual(len(chunks[0]["chunk"]), DEFAULT_CHUNK_WIDTH)

    def test_non_integer_width_falls_back_to_default(self):
        seq = "A" * (DEFAULT_CHUNK_WIDTH + 5)
        chunks = chunk_html(seq, width="not-a-number")
        self.assertEqual(len(chunks), 2)
        self.assertEqual(len(chunks[0]["chunk"]), DEFAULT_CHUNK_WIDTH)

    def test_empty_string_returns_empty_list(self):
        self.assertEqual(chunk_html(""), [])


if __name__ == "__main__":
    unittest.main()
