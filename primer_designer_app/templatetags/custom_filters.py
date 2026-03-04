from django import template

register = template.Library()

@register.filter
def extract_amplicon_info(amplicon_dict, delimiter='|'):
    if len(amplicon_dict["Chrom"].split("|")) == 1:
        return f"{amplicon_dict['Chrom']}: {amplicon_dict['ForPos']} - {amplicon_dict["RevEnd"]}"
    elif len(amplicon_dict["Chrom"].split("|")) > 1:
        transcript_id = amplicon_dict["Chrom"].split("|")[0]
        gene_symbol = amplicon_dict["Chrom"].split("|")[-4]
        
        return f"{gene_symbol} ({transcript_id}): {amplicon_dict['ForPos']} - {amplicon_dict['RevEnd']}"
    else:
        return "Invalid amplicon information"