from django import template
from typing import List, Dict
import re
import logging

LOGGER = logging.getLogger(__name__)

register = template.Library()


@register.filter(name="chunk_seq")
def chunk_seq(seq: str, width: int = 100) -> List[Dict]:
    """
    Split a DNA sequence into chunks of `width` characters.
    Returns list of dicts: {'start': int, 'chunk': str}
    'start' is the 1-based coordinate of the first base of the chunk.
    """
    if seq is None:
        return []
    try:
        width = int(width)
    except Exception:
        width = 100
    chunks = []
    for i in range(0, len(seq), width):
        chunks.append({"start": i + 1, "chunk": seq[i : i + width]})
    return chunks


# HTML-aware chunker: counts visible chars only, preserves tags and re-opens them across lines.
@register.filter(name="chunk_html")
def chunk_html(seq_html: str, width: int = 100) -> List[Dict]:
    """
    Split an HTML string into chunks of `width` visible characters.
    Returns list of dicts: {'start': int, 'chunk': str}
    The returned 'chunk' preserves HTML tags and keeps the markup valid.
    """
    if seq_html is None:
        return []
    try:
        width = int(width)
        if width <= 0:
            width = 100
    except Exception:
        width = 100

    # Tokenize into tags and text
    # matches tags or text: 1. starts with '<' and ends with '>', or 2. any text not containing '<'
    token_re = re.compile(r'(<[^>]+>)|([^<]+)', re.DOTALL)  # matches tags or text
    open_tag_name_re = re.compile(r'^<\s*([A-Za-z0-9:-]+)')
    close_tag_name_re = re.compile(r'^</\s*([A-Za-z0-9:-]+)')
    void_tags = {"br", "img", "hr", "input", "meta", "link", "area", "base", "col", "embed", "param", "source", "track", "wbr"}

    def is_dna_base(ch: str) -> bool:
        return ch.upper() in ("A", "C", "G", "T", "N")

    chunks: List[Dict] = []
    open_tags: List[str] = []          # stack of open tag names
    open_tag_strings: List[str] = []   # corresponding opening tag strings (with attributes)
    current_parts: List[str] = []
    visible_count = 0
    chunk_start = None  # 1-based visible index for current chunk

    def close_all_open_tags() -> str:
        # Return closing tags string in reverse order
        return "".join(f"</{tn}>" for tn in reversed(open_tags))

    def reopen_open_tags() -> List[str]:
        # Return list of opening tag strings to prepend to a new chunk
        return list(open_tag_strings)

    LOGGER.debug(f"Seq HTML: {seq_html}")
    for m in token_re.finditer(seq_html):
        tag = m.group(1)
        text = m.group(2)
        LOGGER.debug(f"Token: tag={tag}, text={text}")
        if tag:
            # handle tag token
            is_closing = tag.startswith("</")
            is_self_closing = tag.endswith("/>") or bool(open_tag_name_re.match(tag) and open_tag_name_re.match(tag).group(1).lower() in void_tags)
            if is_closing:
                # append closing tag
                # try to pop matching open tag from stack
                cname_m = close_tag_name_re.match(tag)
                if cname_m:
                    cname = cname_m.group(1)
                    # pop last matching tag if present
                    for i in range(len(open_tags) - 1, -1, -1):
                        if open_tags[i].lower() == cname.lower():
                            del open_tags[i]
                            del open_tag_strings[i]
                            break
                current_parts.append(tag)
            else:
                # opening or self-closing tag
                if not is_self_closing:
                    nm = open_tag_name_re.match(tag)
                    tname = nm.group(1) if nm else tag
                    open_tags.append(tname)
                    open_tag_strings.append(tag)
                # always append the tag to current parts
                current_parts.append(tag)
            continue

        if text:
            # iterate character-by-character but preserve original chars in output
            pos = 0
            text_len = len(text)
            variant_annotation_ignore = False
            while pos < text_len:
                ch = text[pos]
                if ch == "[":
                    variant_annotation_ignore = True
                elif ch in ["/", ">"]:
                    variant_annotation_ignore = False

                # if starting fresh chunk set chunk_start
                if (visible_count % width) == 0:
                    chunk_start = visible_count + 1
                # append the char to output parts
                current_parts.append(ch)
                # increment visible count only for DNA bases
                if is_dna_base(ch) and not variant_annotation_ignore:
                    visible_count += 1
                pos += 1

                # If we've filled a chunk (visible_count divisible by width) finalize it
                if (visible_count % width == 0) and visible_count > 0 and is_dna_base(ch) and not variant_annotation_ignore:
                    LOGGER.debug(f"Processing char: {ch}, visible_count={visible_count}")
                    # close open tags for this chunk
                    closing = close_all_open_tags()
                    chunk_html_str = "".join(current_parts) + closing
                    chunks.append({"start": chunk_start or 1, "chunk": chunk_html_str})
                    # prepare a new chunk: reopen tags
                    current_parts = reopen_open_tags()[:]  # copy openings
                    # chunk_start becomes next visible index on next iteration
                    chunk_start = None

        LOGGER.debug(f"Current parts: {current_parts}, visible_count: {visible_count}, chunk_start: {chunk_start}")

    # after loop, add remaining content if any visible chars left or tags present
    if current_parts:
        # ensure chunk_start is set for final remainder
        if chunk_start is None:
            chunk_start = visible_count - (visible_count % width) + 1 if visible_count else 1
        closing = close_all_open_tags()
        chunk_html_str = "".join(current_parts) + closing
        chunks.append({"start": chunk_start, "chunk": chunk_html_str})

    LOGGER.debug(f"Final chunks: {chunks}")
    return chunks
