function updateReferenceGenome() {
    const switchButton = document.getElementById("reference-genome-switch");
    const hiddenInput = document.getElementById("reference-genome");
    if (!switchButton || !hiddenInput) return;
    hiddenInput.value = switchButton.checked ? "GRCh38" : "GRCh37";
}

function setProductSizeRange(min, max) {
    const dMin = document.getElementById("dlg_product_size_min");
    const dMax = document.getElementById("dlg_product_size_max");
    if (dMin) dMin.value = min;
    if (dMax) dMax.value = max;
}

/**
 * PCR vs qPCR: sets hidden usecase for the server and applies the usual product
 * size range. Call only when the user explicitly chooses a preset (not on sync/submit).
 */
function applyUsecasePreset(mode) {
    const hidden = document.getElementById("usecase");
    if (hidden) {
        hidden.value = mode === "qPCR" ? "qPCR" : "PCR";
    }
    if (mode === "qPCR") {
        setProductSizeRange(80, 150);
    } else {
        setProductSizeRange(400, 800);
    }
    syncUsecaseButtonHighlight();
}

function syncUsecaseButtonHighlight() {
    const hidden = document.getElementById("usecase");
    if (!hidden) return;
    const mode = hidden.value === "qPCR" ? "qPCR" : "PCR";
    document.querySelectorAll(".primer-usecase-btn").forEach((btn) => {
        const selected = btn.dataset.usecase === mode;
        btn.classList.toggle("is-selected", selected);
        btn.setAttribute("aria-pressed", selected ? "true" : "false");
    });
}

const AMP_CHECK_VALUES = ["none", "genome", "transcriptome"];
const AMP_CHECK_LABELS = ["None", "Genome", "Transcriptome"];

// Update the aria-valuenow and aria-valuetext attributes of the amplicon button
function updateAmpliconCheckAria() {
    const root = document.getElementById("amplicon-button");
    if (!root) return;  // exit if no amplicon button found
    const i = parseInt(root.dataset.position, 10);
    const idx = Number.isNaN(i) ? 0 : i;
    root.setAttribute("aria-valuenow", String(idx));
    root.setAttribute("aria-valuetext", AMP_CHECK_LABELS[idx] || "None");
}

// Set the amplicon check value and update the aria-valuenow
// and aria-valuetext attributes of the amplicon button
function setAmpliconCheck(value) {
    const idx = AMP_CHECK_VALUES.indexOf(value);
    const i = idx >= 0 ? idx : 0;
    const root = document.getElementById("amplicon-button");
    const hidden = document.getElementById("amplicon-check-hidden");
    if (root) root.dataset.position = String(i);
    if (hidden) hidden.value = AMP_CHECK_VALUES[i];
    updateAmpliconCheckAria();
}

// Initialize the amplicon toggle
function initAmpliconToggle() {
    const root = document.getElementById("amplicon-button");
    if (!root) return;
    root.querySelectorAll(".amplicon-opt").forEach((btn) => {
        btn.addEventListener("click", (e) => {
            e.preventDefault();
            const idx = parseInt(btn.dataset.idx, 10);
            if (!Number.isNaN(idx) && idx >= 0 && idx < AMP_CHECK_VALUES.length) {
                setAmpliconCheck(AMP_CHECK_VALUES[idx]);
            }
        });
    });
    root.addEventListener("keydown", (e) => {
        const cur = parseInt(root.dataset.position, 10) || 0;
        if (e.key === "ArrowLeft" && cur > 0) {
            e.preventDefault();
            setAmpliconCheck(AMP_CHECK_VALUES[cur - 1]);
        } else if (e.key === "ArrowRight" && cur < AMP_CHECK_VALUES.length - 1) {
            e.preventDefault();
            setAmpliconCheck(AMP_CHECK_VALUES[cur + 1]);
        }
    });
}

/**
 * Keep hidden POST fields in sync with visible toggles. Browsers may restore
 * checkbox/slider UI from cache while hidden inputs stay at server defaults.
 */
function syncAllTopSettingsFromUi() {
    updateReferenceGenome();
    // Intentionally do not touch usecase or product size here — those come from
    // the primer dialog only; syncing would overwrite user edits before submit.
    const root = document.getElementById("amplicon-button");
    if (root) {
        const pos = parseInt(root.dataset.position, 10);
        if (!Number.isNaN(pos) && pos >= 0 && pos < AMP_CHECK_VALUES.length) {
            setAmpliconCheck(AMP_CHECK_VALUES[pos]);
        }
    }
}

/**
 * Reset primer dialog fields to application defaults for the current use case
 * (PCR: 400–800 bp product; qPCR: 80–150 bp) and Primer3 advanced defaults.
 */
function restorePrimerDialogDefaults() {
    const uc = document.getElementById("usecase");
    const isQ = uc && uc.value === "qPCR";
    const pmin = isQ ? 80 : 400;
    const pmax = isQ ? 150 : 800;
    const setId = (id, val) => {
        const el = document.getElementById(id);
        if (el) el.value = val;
    };
    setId("dlg_tm", "60");
    setId("dlg_gc_content", "50");
    setId("dlg_max_poly_X", "4");
    setId("dlg_product_size_min", String(pmin));
    setId("dlg_product_size_max", String(pmax));

    const dlg = document.getElementById("primer-custom-dialog");
    const p3Defaults = {
        p3_PRIMER_OPT_SIZE: "20",
        p3_PRIMER_MIN_SIZE: "18",
        p3_PRIMER_MAX_SIZE: "22",
        p3_PRIMER_MIN_TM: "",
        p3_PRIMER_MAX_TM: "",
        p3_PRIMER_MIN_GC: "20",
        p3_PRIMER_MAX_GC: "80",
        p3_PRIMER_GC_CLAMP: "1",
        p3_PRIMER_SALT_MONOVALENT: "50",
        p3_PRIMER_DNA_CONC: "50",
        p3_PRIMER_MAX_NS_ACCEPTED: "0",
        p3_PRIMER_MAX_SELF_ANY: "12",
        p3_PRIMER_MAX_SELF_END: "8",
        p3_PRIMER_PAIR_MAX_COMPL_ANY: "12",
        p3_PRIMER_PAIR_MAX_COMPL_END: "8",
        p3_PRIMER_INSIDE_PENALTY: "1",
        p3_PRIMER_INTERNAL_MAX_SELF_END: "8",
        p3_PRIMER_INTERNAL_MAX_POLY_X: "100",
    };
    if (dlg) {
        dlg.querySelectorAll("[name^='p3_']").forEach((el) => {
            const v = p3Defaults[el.name];
            if (v !== undefined) {
                el.value = v;
            }
        });
    }
    syncUsecaseButtonHighlight();
}

function openPrimerCustomDialog() {
    syncUsecaseButtonHighlight();
    const dlg = document.getElementById("primer-custom-dialog");
    if (dlg) dlg.showModal();
}

const identifier = document.querySelector("#Identifier").classList;
const position = document.querySelector("#Genomic_Position").classList;
const sequence = document.querySelector("#Sequence").classList;
const transcriptIDField = document.getElementById("Transcript-ID");
const rel_pos_Field = document.getElementById("Position");
const IdNewBase =  document.getElementById("IDnew_base");
const IdIndelStart = document.getElementById("IdIndelStart");
const IdIndelEnd = document.getElementById("IdIndelEnd");
const IdIndelInst = document.getElementById("IdIndelIns");
const genom_pos_Field = document.getElementById("genom_pos");
const IndelChrom = document.getElementById("IndelChrom");
const IndelStart = document.getElementById("IndelStart");
const IndelEnd = document.getElementById("IndelEnd");
const IndelInst = document.getElementById("IndelIns");
const NewBase =  document.getElementById("new_base");
const sequenceField = document.getElementById("sequence");
const referenceSelectors = Array.from(document.getElementsByName("Reference"));
const cdnaSelector = document.getElementById("cdna");
const cdsSelector = document.getElementById("cds");
const ID_hover_label = document.getElementById("ID-hover-label");
const Genomic_hover_label = document.getElementById("Genomic-hover-label");
const Sequence_hover_label = document.getElementById("Sequence-hover-label");

const reference_genome_switch = document.getElementById("reference-genome-switch");
const switches = [reference_genome_switch];

const IdTranscriptBase = [transcriptIDField, cdnaSelector, cdsSelector];
const IdSNV= [rel_pos_Field, IdNewBase];
const IdSNVList = document.querySelector("#IdSNV").classList;
const IdInDel = [IdIndelStart, IdIndelEnd, IdIndelInst];
const IdInDelList = document.querySelector("#IdInDel").classList;
const idInputFields = [...IdTranscriptBase, ...IdSNV, ...IdInDel];

const SNV = [genom_pos_Field,NewBase];
const SNVList = document.querySelector("#SNV").classList;
const InDel = [IndelChrom, IndelStart, IndelEnd, IndelInst];
const InDelList = document.querySelector("#InDel").classList;
const genomicFields = [...SNV, ...InDel];

function clearAllInputs() {
    const inputs = [...idInputFields, ...genomicFields, sequence, ...switches];
    inputs.forEach(field => {

        field.required = false;
        if (field.type === "radio"){
            field.checked = false;
            console.log(field.checked);
        } else if (field.tagName === "TEXTAREA" || field.tagName === "INPUT") {
            field.value = ""; // Textfelder und Textareas leeren
        } else if (field.tagName === "SELECT") {
            field.selectedIndex = 0; // Dropdown auf die erste Option zurücksetzen
        }
    });
    enableInputField(genomicFields);
    enableInputField(idInputFields);
    enableInputField([sequenceField]);

    document.getElementById("usecase").value = "PCR"; // Standardwert für Use Case
    restorePrimerDialogDefaults();
    document.getElementById("reference-genome").value = "GRCh37"; // Standardwert für Reference Genome
    setAmpliconCheck("none");
    const labels = [ID_hover_label.classList,Genomic_hover_label.classList, Sequence_hover_label.classList]
    labels.forEach(label => label.remove("hover-label-highlight"))

    const containers = [identifier, position, sequence, IdInDelList, IdSNVList, InDelList, SNVList];
    const classesToRemove = ["grey-transparent", "highlight"];
    containers.forEach(el => {
        classesToRemove.forEach(cls => el.remove(cls));
    });
}

let clickCount = 0;
function loadExampleData() {
    clearAllInputs();
    if (clickCount === 0) {
        document.getElementById("genom_pos").value = "ChrY:2655000";
        document.getElementById("new_base").value = "A";
        handleInputChange('genom_pos');

    }
    else if (clickCount === 1){
        document.getElementById("Transcript-ID").value = "ENST00000340855";
        document.getElementById("IdIndelStart").value = "1240";
        document.getElementById("IdIndelEnd").value = "1244";
        document.getElementById("IdIndelIns").value = "GGTC";// Benutzerdefinierte Parameter ausblenden
        document.getElementById("cds").checked = true;
        setAmpliconCheck("genome");
        handleInputChange('Transcript-ID');
    }
    else{
        document.getElementById("sequence").value = "TGTTTCTGATCTTTGTATGCATATGACACTGTAAATGTGTGTGTGTACATATGGTGGGTGTATATATACATGTGTGTACATAAACACATGTACAATTTATTTTATATCTTTACCCCAGTTTGTGATGTAAAATGTTCCTCTGATATGGATAACCATCAAAACAGAGTCAGAAATCACTGTCCTGGATGTCAAAAGTGTTAAAATGACCTTTCAAGGCTCTTCCATTTGAAAGACTTTGTAATCCCAAGATGTTAAGATTCTCAAGAATTTCTGCTTAATATCAGCTGAACACATCTACCAACTGGTAAGATGTTCCCGCCAAACTTGACGGTATATGTGTGTGTACACACACATATGCCATTTTTCAAAATTAAAAAACAAAAGTGCTTTTAGTTGGAAAGAAAGCATCCAAAAACCACTTTTTTGAGTGCATCAATAATTAGCCTAATTTTTACAATAATGAAACTTTTTATCCACTTGAATTTGGGGAGAAAAAGGCCTGTAAGAAAGCCGGACTGGC[-/ATTGCGCAATGC]AGTTAGCTGAAAGTAGCAAAACGTGAGGAAATGCCAAGATCATGGTATTACAGAATGATAGGCCAGATTTTCCAATACCTTAGCTACTCTCTTGTGTTTTTTGTGTATTTAATTGTTTTGTTTTTTAAGAATGTGAATTTATTTTCTCAGGCCACAAGGGTACCACTGAGTAAGCCGGAACTTAAGTTAATGGAATTGGTATCGGACAAGGCTTCCTTCTCCCACATTTGAGAACTTAAATCAGCCCTTAAGGCCTGAAATTCTTCCGTAAGCACTAATTTAGTAAACACTATCCAGGTGGCCTGAATTGAGTTTCCCAAGTTTCCATACAAGTAAGTAAATCAACATCAGAGATACTGACTGCTAAATACCATATCAGAGGAGGTATATTTCTTGTCTTCCAGGAAGCTTCTACACCCATGTGGTATCTAGCTTTAAGGACTCTAAAAGCGGGACACAGGGAGATGGTGCAAAGGATTTCACCTTACAGCTAATTAGGGTGGACCAGAAACTCAGGTACAAAAAGGAGTCTAGGAAGTCAAGGCTGCCTATGAAGAAGCAGAGTCATGAGCCTTGCTGTAGCTTCTTCTCTACTAACCTCTTGTTTAGATGTGGAAGAAATGCCCAGGGTGCCAACAGGGCAGCCAGCAACTTTCCAGAGCAGGTAAGCAAGCTGCCAGTTCAAGGGGCCAGGATCATACAGTAATCAATGATTGTACTGTGATTTGAATTTCTTTGGAGTTCTGTATCTGCAATTGACCATAAGTACCTTGAAGGCAAGGAATTTGTCATTTCTTCTTTCCATCCCTACTACCTCCACTGGTACGTAGTAGGTGCTAGTAAAAGCAAGCAGGACAAAGTAATGAGGAGACACGTGAAAGCACCTGATTCTGTCAAGCCTGTGTATGATCTTCTCATTATGAATTATGACATCTAGAAGAAGGAGATACCTAAGATAGTCAGAATGGGGATGCGGGTGGCGGTGGGGGATGTGGAG"
        handleInputChange('sequence');
        clickCount = -1;
    }

    clickCount++;
}

function enableInputField(fields){
    fields.forEach(field => {

        field.disabled = false;
        let label = document.querySelector(`label[for="${field.id}"]`);
        if (label) {
            label.style.opacity = "1";  // Label normal sichtbar
            label.style.pointerEvents = "auto";  // Label wieder klickbar
        }
        else{
        field.style.opacity = 1;
        field.style.backgroundColor = "";
        }

    })
}

function disableInputField(fields){
    fields.forEach(field => {
        field.disabled = true;
        if (field.type !== "radio"){
            field.value = "";
        }
        let label = document.querySelector(`label[for="${field.id}"]`);
        if (label) {
            label.style.opacity = 0.3;
            label.style.backgroundColor = "#f0f0f"; // Optional: Label ausgrauen
            label.style.pointerEvents = "none";
            label.style.transition= "all 0.2s ease";
        }

    });
}

function isFieldGroupFilled(fields) {
        return fields.some(field => {
            if (field.type === "radio") {
                return field.checked; // Für Radiobuttons: Prüfen, ob einer ausgewählt ist
            }
            return field.value; // Für andere Felder: Prüfen, ob ein Wert vorhanden ist
        });
}

// Funktion, um die anderen Felder zu deaktivieren, wenn eines ausgefüllt ist
function handleInputChange(inputId) {
    if (isFieldGroupFilled(idInputFields)){
        console.log("ID Input Fields filled");
        disableInputField(genomicFields);
        disableInputField([sequenceField]);
        sequence.add("grey-transparent");
        position.add("grey-transparent");
        IdTranscriptBase.forEach(field => field.required = true);
        if(isFieldGroupFilled(IdSNV)){
            IdInDelList.add("grey-transparent");
            IdInDel.forEach(field => field.required = false);
            IdSNV.forEach(field => field.required = true);
            disableInputField(IdInDel);
        }
        else if(isFieldGroupFilled(IdInDel)){
            IdSNVList.add("grey-transparent");
            IdInDel.forEach(field => field.required = true);
            IdSNV.forEach(field => field.required = false);
            IdIndelInst.required = false;
            disableInputField(IdSNV);
        }
        else{
            enableInputField(idInputFields);
            IdSNV.forEach(field => field.required = false);
            IdInDel.forEach(field => field.required = false);

            IdInDelList.remove("grey-transparent");
            IdSNVList.remove("grey-transparent");
        }
    }
    else if (isFieldGroupFilled(genomicFields)){
        disableInputField(idInputFields);
        disableInputField([sequenceField]);
        identifier.add("grey-transparent");
        sequence.add("grey-transparent");
        if (isFieldGroupFilled(InDel)){
            SNVList.add("grey-transparent");
            InDel.forEach(field => field.required = true);
            SNV.forEach(field => field.required = false);
            IndelInst.required = false;
            disableInputField(SNV);
        }
        else if (isFieldGroupFilled(SNV)){
            InDelList.add("grey-transparent");
            InDel.forEach(field => field.required = false);
            SNV.forEach(field => field.required = true);
            disableInputField(InDel);
        }
    }
    else if (isFieldGroupFilled([sequenceField])){
        disableInputField(idInputFields);
        disableInputField(genomicFields);
        position.add("grey-transparent");
        identifier.add("grey-transparent");
    }
    else{
        [idInputFields,genomicFields,[sequenceField]].forEach(container => enableInputField(container));
        [identifier, position, sequence, IdInDelList, IdSNVList, InDelList, SNVList].forEach(container => container.remove("grey-transparent"));
        IdInDel.forEach(field => field.required = false);
        IdSNV.forEach(field => field.required = false);
        IdTranscriptBase.forEach(field => field.required = false);
        referenceSelectors.forEach(radio => radio.checked = false);
        InDel.forEach(field => field.required = false);
        SNV.forEach(field => field.required = false);
        setAmpliconCheck("none");
    }

    ID_hover_label.classList.toggle("hover-label-highlight", isFieldGroupFilled(idInputFields));
    Genomic_hover_label.classList.toggle("hover-label-highlight", isFieldGroupFilled(genomicFields));
    Sequence_hover_label.classList.toggle("hover-label-highlight", isFieldGroupFilled([sequenceField]));

    identifier.toggle("highlight", isFieldGroupFilled(idInputFields));
    position.toggle("highlight", isFieldGroupFilled(genomicFields));
    sequence.toggle("highlight", isFieldGroupFilled([sequenceField]));
    sequenceField.classList.toggle("resize", isFieldGroupFilled([sequenceField]));
}

document.addEventListener("DOMContentLoaded", function () {
    updateAmpliconCheckAria();
    initAmpliconToggle();
    syncAllTopSettingsFromUi();
    syncUsecaseButtonHighlight();
    const settingsForm = document.getElementById("settings");
    if (settingsForm) {
        settingsForm.addEventListener(
            "submit",
            function () {
                syncAllTopSettingsFromUi();
            },
            true
        );
    }
    window.addEventListener("pageshow", function () {
        syncAllTopSettingsFromUi();
        syncUsecaseButtonHighlight();
    });
    const dlg = document.getElementById("primer-custom-dialog");
    if (dlg) {
        dlg.addEventListener("click", function (e) {
            if (e.target === dlg) dlg.close();
        });
    }
});

document.querySelectorAll("td").forEach(td => {
    td.addEventListener("mouseenter", function () {
        const tooltip = this.querySelector(".tooltip");
        if (tooltip) {
            const rect = tooltip.getBoundingClientRect();
            const windowWidth = window.innerWidth;

            // Reset any previous overflow adjustments
            tooltip.removeAttribute("data-overflow");

            // Check for overflow and adjust position
            if (rect.right > windowWidth) {
                tooltip.setAttribute("data-overflow", "right");
            } else if (rect.left < 0) {
                tooltip.setAttribute("data-overflow", "left");
            }
        }
    });
});
