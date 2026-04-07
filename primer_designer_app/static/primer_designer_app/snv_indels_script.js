function updateReferenceGenome() {
    const switchButton = document.getElementById("reference-genome-switch");
    const hiddenInput = document.getElementById("reference-genome");

    // Aktualisiere den Wert basierend auf dem Zustand des Switch-Buttons
    hiddenInput.value = switchButton.checked ? "GRCh38" : "GRCh37";

    console.log("Selected Reference Genome:", hiddenInput.value); // Debugging
}

function setProductSizeRange(min, max) {
    const hMin = document.getElementById("product_size_min");
    const hMax = document.getElementById("product_size_max");
    const dMin = document.getElementById("dlg_product_size_min");
    const dMax = document.getElementById("dlg_product_size_max");
    if (hMin) hMin.value = min;
    if (hMax) hMax.value = max;
    if (dMin) dMin.value = min;
    if (dMax) dMax.value = max;
}

function updateUsecase() {
    const switchButton = document.getElementById("usecase-switch");
    const hiddenInput = document.getElementById("usecase");

    hiddenInput.value = switchButton.checked ? "qPCR" : "PCR";

    if (switchButton.checked) {
        setProductSizeRange(80, 150);
    } else {
        setProductSizeRange(400, 800);
    }
    console.log("Selected Usecase:", hiddenInput.value);
}

function syncPrimerFieldDisabledState() {
    const custom =
        document.getElementById("primer-settings") &&
        document.getElementById("primer-settings").value === "custom";
    document.querySelectorAll("[data-primer-field='default']").forEach((el) => {
        el.disabled = !!custom;
    });
    document.querySelectorAll("[data-primer-field='custom']").forEach((el) => {
        el.disabled = !custom;
    });
}

function syncDialogFromHidden() {
    const pairs = [
        ["tm", "dlg_tm"],
        ["gc_content", "dlg_gc_content"],
        ["max_poly_X", "dlg_max_poly_X"],
        ["product_size_min", "dlg_product_size_min"],
        ["product_size_max", "dlg_product_size_max"],
    ];
    pairs.forEach(([hid, did]) => {
        const h = document.getElementById(hid);
        const d = document.getElementById(did);
        if (h && d) d.value = h.value;
    });
}

function openPrimerCustomDialog() {
    const ps = document.getElementById("primer-settings");
    if (ps) ps.value = "custom";
    syncDialogFromHidden();
    syncPrimerFieldDisabledState();
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
const IdContextSelectors= Array.from(document.getElementsByName("context"));
const cdnaSelector = document.getElementById("cdna");
const cdsSelector = document.getElementById("cds");
const IdGenomeSelector = document.getElementById("genomic");
const IdTranscriptomeSelector = document.getElementById("transcriptomic");
const ID_hover_label = document.getElementById("ID-hover-label");
const Genomic_hover_label = document.getElementById("Genomic-hover-label");
const Sequence_hover_label = document.getElementById("Sequence-hover-label");

const reference_genome_switch = document.getElementById("reference-genome-switch");
const use_case_switch = document.getElementById("usecase-switch");
const switches = [reference_genome_switch, use_case_switch];

const IdAlwaysRequired = [transcriptIDField, cdnaSelector, cdsSelector, IdGenomeSelector, IdTranscriptomeSelector];
const IdSNV= [rel_pos_Field, IdNewBase];
const IdSNVList = document.querySelector("#IdSNV").classList;
const IdInDel = [IdIndelStart, IdIndelEnd, IdIndelInst];
const IdInDelList = document.querySelector("#IdInDel").classList;
const idInputFields = [...IdAlwaysRequired, ...IdSNV, ...IdInDel];

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

    document.getElementById("primer-settings").value = "default";
    document.getElementById("tm").value = 60;
    document.getElementById("gc_content").value = 50;
    document.getElementById("max_poly_X").value = 4;
    setProductSizeRange(400, 800);
    syncPrimerFieldDisabledState();
    document.getElementById("reference-genome").value = "GRCh37"; // Standardwert für Reference Genome
    document.getElementById("usecase").value = "PCR"; // Standardwert für Use Case
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
        document.getElementById("genomic").checked = true;
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
        IdAlwaysRequired.forEach(field => field.required = true);
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
        IdAlwaysRequired.forEach(field => field.required = false);
        IdContextSelectors.forEach(radio => radio.checked = false);
        referenceSelectors.forEach(radio => radio.checked = false);
        InDel.forEach(field => field.required = false);
        SNV.forEach(field => field.required = false);
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
    syncPrimerFieldDisabledState();
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
