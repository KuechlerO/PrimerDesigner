// Primer parameter / overview utilities were extracted into `primer_params_utils.js`.
// This file now contains only SNV/InDel-specific form behaviour.

// SNV/InDel-specific form logic.
// This file is also included on other pages (e.g. structural variants), where
// these elements do not exist. Keep all DOM lookups nullable so the script
// doesn't throw during initial parse, and make SNV/InDel handlers no-op when
// the form isn't present.
const IS_SNV_INDEL_PAGE = Boolean(document.getElementById("Identifier"));

const identifier = document.querySelector("#Identifier")?.classList;
const position = document.querySelector("#Genomic_Position")?.classList;
const sequence = document.querySelector("#Sequence")?.classList;
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
const switches = reference_genome_switch ? [reference_genome_switch] : [];

const IdTranscriptBase = [transcriptIDField, cdnaSelector, cdsSelector].filter(Boolean);
const IdSNV= [rel_pos_Field, IdNewBase].filter(Boolean);
const IdSNVList = document.querySelector("#IdSNV")?.classList;
const IdInDel = [IdIndelStart, IdIndelEnd, IdIndelInst].filter(Boolean);
const IdInDelList = document.querySelector("#IdInDel")?.classList;
const idInputFields = [...IdTranscriptBase, ...IdSNV, ...IdInDel];

const SNV = [genom_pos_Field,NewBase].filter(Boolean);
const SNVList = document.querySelector("#SNV")?.classList;
const InDel = [IndelChrom, IndelStart, IndelEnd, IndelInst].filter(Boolean);
const InDelList = document.querySelector("#InDel")?.classList;
const genomicFields = [...SNV, ...InDel];

function clearAllInputs() {
    if (!IS_SNV_INDEL_PAGE) return;
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

    const tp = document.getElementById("dlg_target_padding");
    if (tp) tp.value = "50";
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
    if (!IS_SNV_INDEL_PAGE) return;
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
    if (!IS_SNV_INDEL_PAGE) return;
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
    if (!IS_SNV_INDEL_PAGE) return;
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
        if (!IS_SNV_INDEL_PAGE) return false;
        return fields.some(field => {
            if (field.type === "radio") {
                return field.checked; // Für Radiobuttons: Prüfen, ob einer ausgewählt ist
            }
            return field.value; // Für andere Felder: Prüfen, ob ein Wert vorhanden ist
        });
}

// Funktion, um die anderen Felder zu deaktivieren, wenn eines ausgefüllt ist
function handleInputChange(inputId) {
    if (!IS_SNV_INDEL_PAGE) return;
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
