function loadStructuralVariantExample() {
    const chr = document.getElementById("sv_chromosome");
    const start = document.getElementById("sv_start_position");
    const end = document.getElementById("sv_end_position");
    const type = document.getElementById("sv_type");

    if (chr) chr.value = "chr12";
    if (start) start.value = "20000000";
    if (end) end.value = "21000000";
    if (type) type.value = "deletion";
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
}
