// Function to highlight a specific region in the HTML string based on the provided start and end positions
// Takes 0-based start and end positions of the primer region in the sequence 
// (excluding HTML tags) and highlights it in the HTML string
function highlightPrimerRegion(html, start, end) {
    let result = '';
    let baseCount = 0;
    let inHtmlTag = false;
    let highlightOpen = false;
    
    // Logic to handle skipping the reference part of DelIns: [REF/ALT]
    let insideBrackets = false;
    let passedSlash = false;

    for (let i = 0; i < html.length; i++) {
        const char = html[i];

        // 1. Handle HTML Tags (don't count these)
        if (char === '<') {
            inHtmlTag = true;
            if (highlightOpen) result += '</span>'; // Close for valid HTML
            result += char;
            continue;
        }
        if (char === '>' && !insideBrackets) {
            inHtmlTag = false;
            result += char;
            if (highlightOpen) result += '<span class="highlight-primer">'; // Re-open
            continue;
        }
        if (inHtmlTag) {
            result += char;
            continue;
        }

        // 2. Logic for Brackets [REF/ALT]
        if (char === '[') {
            insideBrackets = true;
            passedSlash = false;
        }

        // 3. Logic for the Slash separator
        if (insideBrackets && (char === '/' || char === ">")) {
            passedSlash = true;
            result += char;
            continue;
        }

        if (char === ']') {
            insideBrackets = false;
            passedSlash = false;
        }

        // 4. Decide if we count this character
        // We count if it's a letter AND (we are not in brackets OR we are after the slash)
        const isLetter = char.match(/[A-Za-z]/);
        const shouldCount = isLetter && (!insideBrackets || passedSlash);

        if (shouldCount) {
            if (baseCount === start && !highlightOpen) {
                result += '<span class="highlight-primer">';
                highlightOpen = true;
            }

            result += char;
            baseCount++;

            if (baseCount === end && highlightOpen) {
                result += '</span>';
                highlightOpen = false;
            }
        } else {
            // Just append the char (like the 'GGTC' part or brackets) without incrementing baseCount
            result += char;
        }
    }

    if (highlightOpen) result += '</span>';
    return result;
}

// Function to highlight primers in the sequence
function highlightPrimers(primerF_start, primerF_end, primerR_start, primerR_end) {
    document.querySelectorAll(".seq-chunk").forEach(sc => {
        const originalHTML = sc.getAttribute("data-original");

        // Transform chunk start and end as 0-based indices
        const chunkStart = parseInt(sc.previousElementSibling.textContent, 10) - 1;
        const chunkEnd = chunkStart + 99; // Last row may be shorter, but we can ignore that
        console.log("Chunk:", chunkStart, chunkEnd);

        let newHTML = originalHTML;

        // Highlight forward primer if it overlaps this chunk
        if (primerF_start <= chunkEnd && primerF_end >= chunkStart) {
            const relStart = Math.max(0, primerF_start - chunkStart);
            const relEnd = Math.min(sc.textContent.replace(/[^A-Za-z]/g, '').length, primerF_end - chunkStart + 1);
            console.log("Original fwd: ", primerF_start, primerF_end);
            console.log("Highlight Fwd:", relStart, relEnd);
            newHTML = highlightPrimerRegion(newHTML, relStart, relEnd);
        }

        // Highlight reverse primer if it overlaps this chunk
        if (primerR_start <= chunkEnd && primerR_end >= chunkStart) {
            const relStart = Math.max(0, primerR_start - chunkStart);
            const relEnd = Math.min(sc.textContent.replace(/[^A-Za-z]/g, '').length, primerR_end - chunkStart + 1);
            console.log("Original rev: ", primerR_start, primerR_end);
            console.log("Highlight Rev:", relStart, relEnd);
            newHTML = highlightPrimerRegion(newHTML, relStart, relEnd);
        }

        sc.innerHTML = newHTML;
    });
}