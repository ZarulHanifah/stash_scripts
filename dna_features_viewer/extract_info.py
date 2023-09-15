#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 11:40:14 2023

@author: zarul
"""

from bokeh.plotting import figure, output_file, save

gb_path = "/home/zarul/Zarul/Peatlands/status_nanopore/Methylovirgula/bakta/all_annotations/GCF_019343105.1/GCF_019343105.1.gff3"
gff_records = []

with open(gb_path) as in_handle:
    for rec in GFF.parse(in_handle):
        gff_records.append(rec)
        
region_of_interest = [i for i in gff_records if i.id == contig][0]
region_of_interest.features = [i for i in region_of_interest.features if "Is_circular" not in i.qualifiers]
region_of_interest.features = [i for i in region_of_interest.features if "Is_circular" not in i.qualifiers]
region_of_interest_copy = region_of_interest
for i in region_of_interest.features:
    if i.qualifiers["source"][0] not in ["Infernal", "PILER-CR"]:
        try:
            i.qualifiers["product"] = i.qualifiers["locus_tag"]
            i.qualifiers["gene"] = i.qualifiers["locus_tag"]
        except:
            pass
    
record = BiopythonTranslator().translate_record(region_of_interest)
record = record.crop((b_start_p, b_end_p))

selected_labels = [i.html for i in record.features]

df = []
for i in region_of_interest_copy.features:
    # if i.qualifiers["source"][0] not in ["Infernal", "PILER-CR"]:
    try:
        if i.qualifiers["locus_tag"][0] in selected_labels:
            i.qualifiers["start"] = int(i.location._start)
            i.qualifiers["end"] = int(i.location.end)
            i.qualifiers["strand"] = int(i.location.strand)
            i.qualifiers["type"] = i.type
            df.append(i.qualifiers)
    except:
        print(i.qualifiers)
        
df = pd.DataFrame(df)
df.drop(["phase"], axis = "columns", inplace = True)

for i in ["ID", "Name", "locus_tag", "product", "gene", "source"]:
    df[i] = [i[0] for i in df[i]]


df.set_index("ID", inplace = True)

p = record.plot_with_bokeh()
output_file('plot.html', mode='inline')
save(p)
