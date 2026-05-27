// Highlight regions in chunked sequence views (0-based template coordinates).

function _escapeAttr(value) {
    return String(value)
        .replace(/&/g, "&amp;")
        .replace(/"/g, "&quot;")
        .replace(/</g, "&lt;");
}

function _vcfTooltip(hit) {
    const id = hit.id || "VCF variant";
    const alleles = `${hit.ref || "?"} → ${hit.alt || "?"}`;
    const genomic = hit.chrom && hit.pos ? `chr${hit.chrom}:${hit.pos}` : "";
    return genomic ? `${id} (${genomic}) — ${alleles}` : `${id} — ${alleles}`;
}

function _snpTooltip(hit) {
    const id = hit.id || "variant";
    const maf =
        hit.maf != null && hit.maf !== undefined
            ? `${(Number(hit.maf) * 100).toFixed(2)}%`
            : "n/a";
    const alleles = hit.alleles ? ` — ${hit.alleles}` : "";
    return `${id}: MAF ${maf}${alleles}`;
}

function highlightPrimerRegion(html, start, end, cssClass, title, options) {
    const highlightClass = cssClass || "highlight-primer";
    const titleAttr = title ? ` title="${_escapeAttr(title)}"` : "";
    const countInsertedBases =
        options && options.countInsertedBases === true;
    let result = "";
    let baseCount = 0;
    let inHtmlTag = false;
    let highlightOpen = false;

    let insideBrackets = false;
    let countInsideBracket = false;
    let bracketSeenColon = false;
    let bracketCollectRef = false;

    for (let i = 0; i < html.length; i++) {
        const char = html[i];

        if (char === "<") {
            inHtmlTag = true;
            if (highlightOpen) result += "</span>";
            result += char;
            continue;
        }
        if (char === ">" && inHtmlTag) {
            inHtmlTag = false;
            result += char;
            // Only reopen primer highlight outside variant brackets.
            if (highlightOpen && !insideBrackets) {
                result += `<span class="${highlightClass}"${titleAttr}>`;
            }
            continue;
        }
        if (inHtmlTag) {
            result += char;
            continue;
        }

        if (char === "[") {
            insideBrackets = true;
            countInsideBracket = false;
            bracketSeenColon = false;
            bracketCollectRef = false;
        }

        if (
            insideBrackets &&
            char === "/" &&
            i > 0 &&
            html[i - 1] === "+"
        ) {
            countInsideBracket = countInsertedBases;
        }

        if (insideBrackets && char === ":") {
            // Typed bracket like [SNV:G>A] or [DEL:AG/-]
            bracketSeenColon = true;
            bracketCollectRef = true;
        }
        if (insideBrackets && (char === ">" || char === "/")) {
            // Stop counting reference allele bases at delimiter.
            bracketCollectRef = false;
        }

        if (char === "]") {
            insideBrackets = false;
            countInsideBracket = false;
            bracketSeenColon = false;
            bracketCollectRef = false;
        }

        const isLetter = char.match(/[A-Za-z]/);
        // Count template bases 1:1 even when displayed inside HGVS brackets.
        // - reference allele bases inside the bracket count (after ':' in typed form, or immediately after '[' in untyped [G>A])
        // - inserted bases count only for [+/…] when enabled
        const isRefBaseInBracket =
            insideBrackets &&
            isLetter &&
            !countInsideBracket &&
            ((bracketSeenColon && bracketCollectRef) ||
                (!bracketSeenColon &&
                    i > 0 &&
                    html[i - 1] === "[" &&
                    char.match(/[ACGTNacgtn]/)));

        // WT [-/…]: inserted bases are not on the template. MUT [+/…]: count them.
        const shouldCount =
            isLetter && (!insideBrackets || countInsideBracket || isRefBaseInBracket);

        if (shouldCount) {
            if (baseCount === start && !highlightOpen) {
                result += `<span class="${highlightClass}"${titleAttr}>`;
                highlightOpen = true;
            }

            result += char;
            baseCount++;

            if (baseCount === end && highlightOpen) {
                result += "</span>";
                highlightOpen = false;
            }
        } else {
            result += char;
        }
    }

    if (highlightOpen) result += "</span>";
    return result;
}

function _chunkTemplateLetterLength(html, options) {
    const countInsertedBases =
        options && options.countInsertedBases === true;
    let count = 0;
    let inHtmlTag = false;
    let insideBrackets = false;
    let countInsideBracket = false;
    let bracketSeenColon = false;
    let bracketCollectRef = false;

    for (let i = 0; i < html.length; i++) {
        const char = html[i];
        if (char === "<") {
            inHtmlTag = true;
            continue;
        }
        if (char === ">") {
            inHtmlTag = false;
            continue;
        }
        if (inHtmlTag) {
            continue;
        }
        if (char === "[") {
            insideBrackets = true;
            countInsideBracket = false;
            bracketSeenColon = false;
            bracketCollectRef = false;
        }
        if (
            insideBrackets &&
            char === "/" &&
            i > 0 &&
            html[i - 1] === "+"
        ) {
            countInsideBracket = countInsertedBases;
        }
        if (insideBrackets && char === ":") {
            bracketSeenColon = true;
            bracketCollectRef = true;
        }
        if (insideBrackets && (char === ">" || char === "/")) {
            bracketCollectRef = false;
        }
        if (char === "]") {
            insideBrackets = false;
            countInsideBracket = false;
            bracketSeenColon = false;
            bracketCollectRef = false;
        }
        const isLetter = char.match(/[A-Za-z]/);
        const isRefBaseInBracket =
            insideBrackets &&
            isLetter &&
            !countInsideBracket &&
            ((bracketSeenColon && bracketCollectRef) ||
                (!bracketSeenColon &&
                    i > 0 &&
                    html[i - 1] === "[" &&
                    char.match(/[ACGTNacgtn]/)));
        if (isLetter && (!insideBrackets || countInsideBracket || isRefBaseInBracket)) {
            count++;
        }
    }
    return count;
}

function _isDnaLetter(ch) {
    return ch && ch.match(/[ACGTNacgtn]/);
}

/** How many template bases a HGVS bracket replaces (SNV=1, INS=0, DEL=len(ref), …). */
function _templateBasesConsumedByBracket(innerPlain) {
    let body = innerPlain;
    const colon = body.indexOf(":");
    if (colon >= 0) {
        body = body.slice(colon + 1);
    }
    if (body.includes(">")) {
        return 1;
    }
    if (body.startsWith("-/")) {
        return 0;
    }
    if (body.endsWith("/-")) {
        return body.slice(0, body.indexOf("/-")).length;
    }
    const slash = body.indexOf("/");
    if (slash >= 0) {
        return body.slice(0, slash).length;
    }
    return 1;
}

/**
 * Apply primer highlights on annotated HTML using plain-template coordinates.
 * Walks plain and annotated in parallel so HGVS brackets cannot desync indices.
 */
function applyPrimersByPlainCoords(annotatedHtml, plainChunk, ranges) {
    let result = "";
    let plainIdx = 0;
    let highlightOpen = false;
    let activeTitle = "";
    const plainLen = plainChunk.length;

    const activeRange = (idx) => {
        for (const r of ranges) {
            if (idx >= r.start && idx < r.end) {
                return r;
            }
        }
        return null;
    };

    const openHighlight = (range) => {
        if (!range) {
            return;
        }
        const titleAttr = range.title
            ? ` title="${_escapeAttr(range.title)}"`
            : "";
        if (!highlightOpen || activeTitle !== range.title) {
            if (highlightOpen) {
                result += "</span>";
            }
            result += `<span class="${range.cssClass || "highlight-primer"}"${titleAttr}>`;
            highlightOpen = true;
            activeTitle = range.title || "";
        }
    };
    const closeHighlight = () => {
        if (highlightOpen) {
            result += "</span>";
            highlightOpen = false;
            activeTitle = "";
        }
    };

    const emitTemplateBase = (templateIdx, ch) => {
        const range = activeRange(templateIdx);
        if (range) {
            openHighlight(range);
        } else {
            closeHighlight();
        }
        result += ch;
    };

    let ai = 0;
    while (ai < annotatedHtml.length) {
        if (annotatedHtml[ai] === "<") {
            closeHighlight();
            const tagStart = ai;
            while (ai < annotatedHtml.length && annotatedHtml[ai] !== ">") {
                ai++;
            }
            if (ai < annotatedHtml.length) {
                result += annotatedHtml.slice(tagStart, ai + 1);
                ai++;
            } else {
                result += annotatedHtml.slice(tagStart);
                break;
            }
            const range = activeRange(plainIdx);
            if (range) {
                openHighlight(range);
            }
            continue;
        }

        if (annotatedHtml[ai] === "[") {
            closeHighlight();
            const bracketPlainStart = plainIdx;
            let bracketEnd = ai;
            while (
                bracketEnd < annotatedHtml.length &&
                annotatedHtml[bracketEnd] !== "]"
            ) {
                bracketEnd++;
            }
            if (bracketEnd >= annotatedHtml.length) {
                result += annotatedHtml.slice(ai);
                break;
            }

            let innerPlain = "";
            let bi = ai + 1;
            while (bi < bracketEnd) {
                if (annotatedHtml[bi] === "<") {
                    bi++;
                    while (bi < bracketEnd && annotatedHtml[bi] !== ">") {
                        bi++;
                    }
                    if (bi < bracketEnd) {
                        bi++;
                    }
                    continue;
                }
                innerPlain += annotatedHtml[bi];
                bi++;
            }
            const isSnvBracket = innerPlain.includes(">");

            // For non-SNV brackets, the number of template bases consumed depends on
            // which allele template we're rendering (WT has ref bases, MUT has alt bases).
            let consumed = _templateBasesConsumedByBracket(innerPlain);
            if (!isSnvBracket) {
                let insBody = innerPlain;
                const insColon = insBody.indexOf(":");
                if (insColon >= 0) {
                    insBody = insBody.slice(insColon + 1);
                }
                if (insBody.startsWith("-/")) {
                    const insSeq = insBody.slice(2);
                    const plainFromHere = plainChunk.slice(bracketPlainStart);
                    if (insSeq && plainFromHere.startsWith(insSeq)) {
                        consumed = insSeq.length;
                    } else {
                        consumed = 0;
                    }
                }
            }
            let bracketBody = innerPlain;
            const colonIdx = bracketBody.indexOf(":");
            if (colonIdx >= 0) {
                bracketBody = bracketBody.slice(colonIdx + 1);
            }
            let refPart = "";
            let altPart = "";
            let allelePart = "ref"; // which side matches the current plain template
            if (!isSnvBracket) {
                const slashIdx = bracketBody.indexOf("/");
                if (slashIdx >= 0) {
                    refPart = bracketBody.slice(0, slashIdx);
                    altPart = bracketBody.slice(slashIdx + 1);
                    const plainFromHere = plainChunk.slice(bracketPlainStart);
                    if (refPart && plainFromHere.startsWith(refPart)) {
                        allelePart = "ref";
                        consumed = refPart.length;
                    } else if (altPart && plainFromHere.startsWith(altPart)) {
                        allelePart = "alt";
                        consumed = altPart.length;
                    }
                }
            }

            let refAlleleSeen = false;
            let afterColon = false;
            let bracketPlainCursor = bracketPlainStart;
            let snvPhase = "prefix";
            for (let k = ai; k <= bracketEnd; k++) {
                const ch = annotatedHtml[k];
                if (ch === "<") {
                    closeHighlight();
                    const ts = k;
                    let te = k;
                    while (te < annotatedHtml.length && annotatedHtml[te] !== ">") {
                        te++;
                    }
                    result += annotatedHtml.slice(ts, te + 1);
                    k = te;
                    continue;
                }
                if (ch === ":") {
                    afterColon = true;
                    if (isSnvBracket) {
                        snvPhase = "ref";
                    }
                    closeHighlight();
                    result += ch;
                    continue;
                }
                if (ch === "/" && !isSnvBracket) {
                    closeHighlight();
                    result += ch;
                    continue;
                }
                if (ch === "]") {
                    closeHighlight();
                    result += ch;
                    continue;
                }
                if (ch === ">") {
                    if (isSnvBracket) {
                        snvPhase = "alt";
                    }
                    closeHighlight();
                    result += ch;
                    continue;
                }
                if (_isDnaLetter(ch)) {
                    let templateIdx = -1;
                    if (isSnvBracket) {
                        // Only ref/alt letters (not "N" in "SNV") may map to template index.
                        if (
                            (snvPhase === "ref" || snvPhase === "alt") &&
                            bracketPlainStart < plainLen &&
                            ch.toUpperCase() ===
                                plainChunk[bracketPlainStart].toUpperCase()
                        ) {
                            templateIdx = bracketPlainStart;
                        }
                    } else {
                        // For DEL/DELINS notation: map every base from the allele-present side
                        // (WT uses the ref part before '/', MUT uses the alt part after '/').
                        let phase = "prefix";
                        if (afterColon) {
                            phase = "ref";
                        }
                        // Look back to see if we already crossed the slash in this bracket.
                        for (let back = k - 1; back >= ai; back--) {
                            if (annotatedHtml[back] === "]") break;
                            if (annotatedHtml[back] === "/") {
                                phase = "alt";
                                break;
                            }
                            if (annotatedHtml[back] === ":") {
                                break;
                            }
                        }

                        const eligible =
                            (allelePart === "ref" && phase === "ref") ||
                            (allelePart === "alt" && phase === "alt");
                        if (eligible && bracketPlainCursor < bracketPlainStart + consumed) {
                            templateIdx = bracketPlainCursor;
                            bracketPlainCursor += 1;
                        } else {
                            // Only fall back to the legacy "first ref base" behavior if this bracket
                            // is NOT a ref/alt bracket (i.e. no slash present).
                            if (!refPart && !altPart) {
                                const isRefAllele =
                                    !refAlleleSeen &&
                                    (afterColon ||
                                        (k > ai && annotatedHtml[k - 1] === "["));
                                if (isRefAllele) {
                                    templateIdx = bracketPlainCursor;
                                    refAlleleSeen = true;
                                    bracketPlainCursor += 1;
                                }
                            }
                        }
                    }
                    if (templateIdx >= 0) {
                        emitTemplateBase(templateIdx, ch);
                    } else {
                        closeHighlight();
                        result += ch;
                    }
                    continue;
                }
                closeHighlight();
                result += ch;
            }

            plainIdx = bracketPlainStart + consumed;
            ai = bracketEnd + 1;
            if (!activeRange(plainIdx)) {
                closeHighlight();
            }
            continue;
        }

        if (plainIdx < plainLen && _isDnaLetter(plainChunk[plainIdx])) {
            const ch = annotatedHtml[ai];
            if (_isDnaLetter(ch) && ch.toUpperCase() === plainChunk[plainIdx].toUpperCase()) {
                emitTemplateBase(plainIdx, ch);
                plainIdx++;
                ai++;
                if (!activeRange(plainIdx)) {
                    closeHighlight();
                }
                continue;
            }
        }

        closeHighlight();
        result += annotatedHtml[ai];
        ai++;
    }

    closeHighlight();
    return result;
}

/** Restore variant highlight after primer layers so red stays visible under teal. */
function _reapplyMutationHighlights(html, originalHTML) {
    const mutRe = /<span class=['"]highlight-mutation['"][^>]*>([\s\S]*?)<\/span>/i;
    const mutMatch = originalHTML.match(mutRe);
    if (!mutMatch) {
        return html;
    }
    const variantText = mutMatch[1].replace(/<[^>]+>/g, "");
    if (!variantText) {
        return html;
    }
    // Guard: if we ever captured a whole line, don't repaint it.
    if (variantText.length > 60) {
        return html;
    }

    let result = "";
    let vi = 0;
    let inHtmlTag = false;
    let mutationOpen = false;

    for (let i = 0; i < html.length; i++) {
        const char = html[i];
        if (char === "<") {
            inHtmlTag = true;
            if (mutationOpen) {
                result += "</span>";
                mutationOpen = false;
            }
            result += char;
            continue;
        }
        if (char === ">" && inHtmlTag) {
            inHtmlTag = false;
            result += char;
            continue;
        }
        if (inHtmlTag) {
            result += char;
            continue;
        }

        if (vi < variantText.length && char === variantText[vi]) {
            if (!mutationOpen) {
                result += "<span class='highlight-mutation'>";
                mutationOpen = true;
            }
            result += char;
            vi++;
            if (vi >= variantText.length && mutationOpen) {
                result += "</span>";
                mutationOpen = false;
            }
        } else {
            if (mutationOpen) {
                result += "</span>";
                mutationOpen = false;
            }
            result += char;
        }
    }
    if (mutationOpen) {
        result += "</span>";
    }
    return result;
}

function _isHitBindingConflict(
    hit,
    primerF_start,
    primerF_end,
    primerR_start,
    primerR_end
) {
    const ts = hit.template_start;
    const te = hit.template_end;
    const inForward = ts <= primerF_end && te >= primerF_start;
    const inReverse = ts <= primerR_end && te >= primerR_start;
    return inForward || inReverse;
}

/**
 * Rebuild sequence chunk HTML from data-original, applying layers in order:
 * VCF spikes (lavender) → common SNPs (amber) → binding-site SNPs (orange) → primers (teal).
 */
function _applyVcfHighlights(html, chunkStart, chunkEnd) {
    const hits =
        typeof VCF_REGION_HITS !== "undefined" && VCF_REGION_HITS.length
            ? VCF_REGION_HITS
            : [];
    for (const hit of hits) {
        const ts = hit.template_start;
        const te = hit.template_end;
        if (ts > chunkEnd || te < chunkStart) {
            continue;
        }
        const relStart = Math.max(0, ts - chunkStart);
        const relEnd = Math.min(
            chunkEnd - chunkStart + 1,
            te - chunkStart + 1
        );
        html = highlightPrimerRegion(
            html,
            relStart,
            relEnd,
            "highlight-vcf",
            _vcfTooltip(hit)
        );
    }
    return html;
}

function refreshSequenceHighlights(
    primerF_start,
    primerF_end,
    primerR_start,
    primerR_end
) {
    const root = arguments.length >= 5 && arguments[4] ? arguments[4] : document;
    const countInsertedBases =
        root.getAttribute &&
        root.getAttribute("data-count-bracket-insertion") === "1";
    const isAlleleSpecific =
        root.getAttribute && root.getAttribute("data-allele-specific") === "1";
    const primerHighlightOpts = { countInsertedBases };
    const hits =
        typeof SNP_REGION_HITS !== "undefined" && SNP_REGION_HITS.length
            ? SNP_REGION_HITS
            : [];

    root.querySelectorAll(".seq-chunk").forEach((sc) => {
        const originalHTML = sc.getAttribute("data-original");
        const annotatedHTML = sc.getAttribute("data-annotated") || originalHTML;
        const chunkStart = parseInt(sc.previousElementSibling.textContent, 10) - 1;
        const chunkLen = isAlleleSpecific
            ? originalHTML.length
            : _chunkTemplateLetterLength(originalHTML, primerHighlightOpts);
        const chunkEnd = chunkStart + Math.max(0, chunkLen - 1);

        let html = isAlleleSpecific ? annotatedHTML : originalHTML;
        if (!isAlleleSpecific) {
            html = _applyVcfHighlights(html, chunkStart, chunkEnd);
        }

        const primerLayers = [
            ["F", primerF_start, primerF_end, "highlight-primer"],
            ["R", primerR_start, primerR_end, "highlight-primer"],
        ];
        const primerRanges = [];
        for (const [pLabel, pStart, pEnd, cssClass] of primerLayers) {
            if (pStart <= chunkEnd && pEnd >= chunkStart) {
                primerRanges.push({
                    start: Math.max(0, pStart - chunkStart),
                    end: Math.min(
                        isAlleleSpecific
                            ? originalHTML.length
                            : _chunkTemplateLetterLength(html, primerHighlightOpts),
                        pEnd - chunkStart + 1
                    ),
                    cssClass,
                    title: `Primer ${pLabel}`,
                });
            }
        }
        if (isAlleleSpecific && primerRanges.length) {
            html = applyPrimersByPlainCoords(
                annotatedHTML,
                originalHTML,
                primerRanges
            );
        } else if (!isAlleleSpecific) {
            for (const range of primerRanges) {
                html = highlightPrimerRegion(
                    html,
                    range.start,
                    range.end,
                    range.cssClass,
                    range.title,
                    primerHighlightOpts
                );
            }
        }

        // Reapply variant highlight from annotated HTML (SNV/Indel views only).
        if (!isAlleleSpecific) {
            html = _reapplyMutationHighlights(html, annotatedHTML);
        }

        const nonConflictHits = [];
        const conflictHits = [];
        for (const hit of hits) {
            const ts = hit.template_start;
            const te = hit.template_end;
            if (ts > chunkEnd || te < chunkStart) {
                continue;
            }
            if (
                _isHitBindingConflict(
                    hit,
                    primerF_start,
                    primerF_end,
                    primerR_start,
                    primerR_end
                )
            ) {
                conflictHits.push(hit);
            } else {
                nonConflictHits.push(hit);
            }
        }

        for (const hit of nonConflictHits) {
            const relStart = Math.max(0, hit.template_start - chunkStart);
            const relEnd = Math.min(
                chunkEnd - chunkStart + 1,
                hit.template_end - chunkStart + 1
            );
            html = highlightPrimerRegion(
                html,
                relStart,
                relEnd,
                "highlight-snp",
                _snpTooltip(hit)
            );
        }

        for (const hit of conflictHits) {
            const relStart = Math.max(0, hit.template_start - chunkStart);
            const relEnd = Math.min(
                chunkEnd - chunkStart + 1,
                hit.template_end - chunkStart + 1
            );
            html = highlightPrimerRegion(
                html,
                relStart,
                relEnd,
                "highlight-snp-conflict",
                _snpTooltip(hit)
            );
        }

        sc.innerHTML = html;
    });
}

function highlightAllSnpsInRegion(root) {
    const scope = root || document;
    scope.querySelectorAll(".seq-chunk").forEach((sc) => {
        const originalHTML = sc.getAttribute("data-original");
        const chunkStart = parseInt(sc.previousElementSibling.textContent, 10) - 1;
        const countInsertedBases =
            scope.getAttribute &&
            scope.getAttribute("data-count-bracket-insertion") === "1";
        const primerHighlightOpts = { countInsertedBases };
        const chunkLen = _chunkTemplateLetterLength(originalHTML, primerHighlightOpts);
        const chunkEnd = chunkStart + Math.max(0, chunkLen - 1);
        let html = originalHTML;
        html = _applyVcfHighlights(html, chunkStart, chunkEnd);
        const hits =
            typeof SNP_REGION_HITS !== "undefined" ? SNP_REGION_HITS : [];
        for (const hit of hits) {
            const ts = hit.template_start;
            const te = hit.template_end;
            if (ts > chunkEnd || te < chunkStart) {
                continue;
            }
            const relStart = Math.max(0, ts - chunkStart);
            const relEnd = Math.min(
                chunkEnd - chunkStart + 1,
                te - chunkStart + 1
            );
            html = highlightPrimerRegion(
                html,
                relStart,
                relEnd,
                "highlight-snp",
                _snpTooltip(hit)
            );
        }
        sc.innerHTML = html;
    });
}

function _sequenceHighlightRoot() {
    return (
        document.getElementById("sequence-container-wt") ||
        document.getElementById("sequence-container-mut") ||
        document.getElementById("sequence-container") ||
        document
    );
}

function highlightSnpsForPair(primerF_start, primerF_end, primerR_start, primerR_end) {
    refreshSequenceHighlights(
        primerF_start,
        primerF_end,
        primerR_start,
        primerR_end,
        _sequenceHighlightRoot()
    );
}

function highlightPrimers(primerF_start, primerF_end, primerR_start, primerR_end) {
    refreshSequenceHighlights(
        primerF_start,
        primerF_end,
        primerR_start,
        primerR_end,
        _sequenceHighlightRoot()
    );
}
