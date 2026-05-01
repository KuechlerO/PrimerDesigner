function loadStructuralVariantExample() {
    const chr = document.getElementById("sv_chromosome");
    const start = document.getElementById("sv_start_position");
    const end = document.getElementById("sv_end_position");
    const type = document.getElementById("sv_type");

    if (chr) chr.value = "chr12";
    if (start) start.value = "20000000";
    if (end) end.value = "21000000";
    if (type) type.value = "deletion";
    updateSvSketch();
}

function clearStructuralVariantInputs() {
    const chr = document.getElementById("sv_chromosome");
    const start = document.getElementById("sv_start_position");
    const end = document.getElementById("sv_end_position");
    const type = document.getElementById("sv_type");

    if (chr) chr.value = "";
    if (start) start.value = "";
    if (end) end.value = "";
    if (type) type.value = "deletion";
    updateSvSketch();
}

function parsePos(value) {
    if (value == null) return null;
    const cleaned = String(value).trim().replaceAll(",", "").replaceAll("_", "");
    if (!cleaned) return null;
    const n = Number(cleaned);
    return Number.isFinite(n) ? n : null;
}

function formatInt(n) {
    try {
        return new Intl.NumberFormat(undefined, { maximumFractionDigits: 0 }).format(n);
    } catch {
        return String(n);
    }
}

function svgEl(tag, attrs = {}) {
    const el = document.createElementNS("http://www.w3.org/2000/svg", tag);
    for (const [k, v] of Object.entries(attrs)) {
        if (v === undefined || v === null) continue;
        el.setAttribute(k, String(v));
    }
    return el;
}

function clearSvg(svg) {
    while (svg.firstChild) svg.removeChild(svg.firstChild);
}

function updateSvSketch() {
    const svg = document.getElementById("sv_sketch_svg");
    if (!svg) return;

    const chr = document.getElementById("sv_chromosome")?.value?.trim() ?? "";
    const start = parsePos(document.getElementById("sv_start_position")?.value);
    const end = parsePos(document.getElementById("sv_end_position")?.value);

    clearSvg(svg);

    // Layout constants (match viewBox="0 0 520 170")
    const W = 520;
    const H = 170;
    const padX = 34;
    const axisY = 92;
    const axisX0 = padX;
    const axisX1 = W - padX;

    const roiValid = start != null && end != null && end > start;
    const roiLen = roiValid ? (end - start) : null;

    // Choose a design window around ROI for the sketch
    const flank = roiValid ? Math.max(Math.round(roiLen * 0.35), 1000) : 5000;
    const winStart = roiValid ? (start - flank) : 0;
    const winEnd = roiValid ? (end + flank) : 1;
    const scale = (axisX1 - axisX0) / (winEnd - winStart);
    const xFor = (pos) => axisX0 + (pos - winStart) * scale;

    // Background
    svg.appendChild(svgEl("rect", {
        x: 0, y: 0, width: W, height: H,
        rx: 10, ry: 10,
        fill: "rgba(255,255,255,0.0)"
    }));

    // Axis line
    svg.appendChild(svgEl("line", {
        x1: axisX0, y1: axisY, x2: axisX1, y2: axisY,
        stroke: "rgb(87, 163, 176)",
        "stroke-width": 4,
        "stroke-linecap": "round",
        opacity: 0.85
    }));

    // Axis end ticks
    for (const x of [axisX0, axisX1]) {
        svg.appendChild(svgEl("line", {
            x1: x, y1: axisY - 10, x2: x, y2: axisY + 10,
            stroke: "rgba(39,52,56,0.35)",
            "stroke-width": 2
        }));
    }

    // Labels
    const titleText = roiValid
        ? `${chr || "chr?"}: ${formatInt(start)} – ${formatInt(end)}`
        : (chr ? `${chr}: (enter start/end)` : "Enter chromosome/start/end");

    svg.appendChild(svgEl("text", {
        x: W / 2,
        y: 26,
        "text-anchor": "middle",
        "font-family": "system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif",
        "font-size": 16,
        "font-weight": 700,
        fill: "rgb(45, 133, 120)"
    })).textContent = titleText;

    // ROI highlight
    if (roiValid) {
        const x0 = xFor(start);
        const x1 = xFor(end);

        svg.appendChild(svgEl("rect", {
            x: x0,
            y: axisY - 18,
            width: Math.max(2, x1 - x0),
            height: 36,
            rx: 8,
            fill: "rgba(82, 175, 162, 0.22)",
            stroke: "rgba(82, 175, 162, 0.8)",
            "stroke-width": 2
        }));

        svg.appendChild(svgEl("text", {
            x: (x0 + x1) / 2,
            y: H - 8,
            "text-anchor": "middle",
            "font-family": "system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif",
            "font-size": 12,
            fill: "rgba(39,52,56,0.75)"
        })).textContent = `ROI (${formatInt(end - start)} bp)`;
    }

    function drawPrimerPair(xMid, label, color, yOffset = 0, labelPlacement = "above") {
        const y = axisY + yOffset;
        const halfSpan = 34;
        const xL = xMid - halfSpan;
        const xR = xMid + halfSpan;

        const arrowLen = 10;
        const arrowHalfHeight = 5;

        // left primer
        svg.appendChild(svgEl("polygon", {
            points: `${xL},${y} ${xL + arrowLen},${y - arrowHalfHeight} ${xL + arrowLen},${y + arrowHalfHeight}`,
            fill: color,
            opacity: 0.95
        }));
        // right primer
        svg.appendChild(svgEl("polygon", {
            points: `${xR},${y} ${xR - arrowLen},${y - arrowHalfHeight} ${xR - arrowLen},${y + arrowHalfHeight}`,
            fill: color,
            opacity: 0.95
        }));
        // amplicon line
        svg.appendChild(svgEl("line", {
            x1: xL + arrowLen,
            y1: y,
            x2: xR - arrowLen,
            y2: y,
            stroke: color,
            "stroke-width": 3,
            "stroke-linecap": "round",
            opacity: 0.9
        }));

        const labelY = labelPlacement === "below" ? (y + 18) : (y - 14);
        svg.appendChild(svgEl("text", {
            x: xMid,
            y: labelY,
            "text-anchor": "middle",
            "font-family": "system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif",
            "font-size": 11,
            "font-weight": 700,
            fill: "rgba(39,52,56,0.72)"
        })).textContent = label;
    }

    // Primer pairs: 1 upstream, 2 inside ROI, 1 downstream
    const cUp = "rgb(87, 163, 176)";
    const cIn = "rgb(45, 133, 120)";
    const cDown = "rgb(153, 199, 207)";

    if (roiValid) {
        const xUp = xFor(start - Math.max(Math.round(roiLen * 0.18), 800));
        const xIn1 = xFor(start + Math.round(roiLen * 0.30));
        const xIn2 = xFor(start + Math.round(roiLen * 0.70));
        const xDown = xFor(end + Math.max(Math.round(roiLen * 0.18), 800));

        drawPrimerPair(xUp, "Upstream", cUp, -26, "above");
        // More distance from the main axis for Internal 1
        drawPrimerPair(xIn1, "Internal 1", cIn, -18, "above");
        // Place Internal 2 on ROI border (not inside the ROI box)
        drawPrimerPair(xIn2, "Internal 2", cIn, 18, "below");
        drawPrimerPair(xDown, "Downstream", cDown, 34, "below");

        // Start/end ticks
        for (const [pos, tLabel] of [[start, "start"], [end, "end"]]) {
            const x = xFor(pos);
            svg.appendChild(svgEl("line", {
                x1: x, y1: axisY - 26, x2: x, y2: axisY + 26,
                stroke: "rgba(39,52,56,0.35)",
                "stroke-width": 2
            }));
            svg.appendChild(svgEl("text", {
                x,
                y: axisY - 32,
                "text-anchor": "middle",
                "font-family": "system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif",
                "font-size": 11,
                fill: "rgba(39,52,56,0.65)"
            })).textContent = tLabel;
        }

    } else {
        // Placeholder “expected” layout (static positions)
        drawPrimerPair(axisX0 + (axisX1 - axisX0) * 0.2, "Upstream", cUp, -18, "above");
        drawPrimerPair(axisX0 + (axisX1 - axisX0) * 0.45, "Internal 1", cIn, -6, "above");
        drawPrimerPair(axisX0 + (axisX1 - axisX0) * 0.55, "Internal 2", cIn, 24, "below");
        drawPrimerPair(axisX0 + (axisX1 - axisX0) * 0.8, "Downstream", cDown, 44, "below");
    }
}

document.addEventListener("DOMContentLoaded", () => {
    // Structural variant page default: qPCR preset (only if user hasn't changed defaults)
    // This relies on applyAssayPreset() from primer_params_utils.js.
    const padEl = document.getElementById("dlg_target_padding");
    const pminEl = document.getElementById("dlg_product_size_min");
    const pmaxEl = document.getElementById("dlg_product_size_max");
    const pad = parseInt(padEl?.value, 10);
    const pmin = parseInt(pminEl?.value, 10);
    const pmax = parseInt(pmaxEl?.value, 10);
    const looksLikePcrDefaults = pad === 50 && pmin === 400 && pmax === 800;
    if (looksLikePcrDefaults && typeof applyAssayPreset === "function") {
        applyAssayPreset("qPCR");
    } else if (typeof syncAssayPresetButtonHighlight === "function") {
        syncAssayPresetButtonHighlight();
    }

    const chr = document.getElementById("sv_chromosome");
    const start = document.getElementById("sv_start_position");
    const end = document.getElementById("sv_end_position");
    [chr, start, end].forEach((el) => {
        if (!el) return;
        el.addEventListener("input", updateSvSketch);
        el.addEventListener("change", updateSvSketch);
    });
    updateSvSketch();
});
