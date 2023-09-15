#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 16:57:58 2023

@author: zarul
"""

from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt
from BCBio import GFF

# plt.style.use("ggplot")

def gc_content(s):
    return 100.0 * len([c for c in s if c in "GC"]) / len(s)

def gen_local_gc_content(seq, window_size):
    yy = [
        gc_content(seq[i - int(window_size/2) : i + int(window_size/2)])
        for i in range(int(window_size/2), len(seq) - int(window_size/2))
    ]
    return yy

final_record = {}

# Sulfur-related genes bound by IS5 transposase, position 54150 until 102423
# title = "Sulfur genes"
# contig, b_start, b_end = "contig_1", 54150, 102423

# Methane-related genes are bound by IS66 transposase, position 3010377 until 3043313
title = "Methane genes"
contig, b_start, b_end = "contig_1", 3010377, 3043313

padding = int(round((b_end - b_start)/2*0.1))

b_start_p, b_end_p = b_start - padding, b_end + padding

fig, axs = plt.subplots(nrows = 3, figsize = (10, 5), dpi = 300)

# ax 0
gb_path = "/home/zarul/Zarul/Peatlands/status_nanopore/Methylovirgula/bakta/all_annotations/GCF_019343105.1/GCF_019343105.1.gff3"
gff_records = []

with open(gb_path) as in_handle:
    for rec in GFF.parse(in_handle):
        gff_records.append(rec)

region_of_interest = [i for i in gff_records if i.id == contig][0]
seq = region_of_interest._seq
region_of_interest.features = [i for i in region_of_interest.features if "Is_circular" not in i.qualifiers] # remove contig features

record = BiopythonTranslator().translate_record(region_of_interest)
record = record.crop((b_start_p, b_end_p))
final_record["bakta"] = record

for feature in record.features:
    feature.label = None

record.plot(ax = axs[0], x_lim=(b_start_p, b_end_p))
axs[0].annotate("Prodigal-Bakta", (b_start_p - padding, 0), rotation = 90, annotation_clip=False)

# ax 1
gb_path = "/home/zarul/Zarul/Peatlands/status_nanopore/Methylovirgula/Gwak_isescane/GCF_019343105.1/GCF_019343105.1.fna.gff"
gff_records = []

with open(gb_path) as in_handle:
    for rec in GFF.parse(in_handle):
        gff_records.append(rec)

region_of_interest = [i for i in gff_records if i.id == contig][0]
region_of_interest._seq = seq
region_of_interest.features = [i for i in region_of_interest.features if "Is_circular" not in i.qualifiers] # remove contig features

for i in region_of_interest.features:
    try:
        i.qualifiers["source"] = i.qualifiers["family"]
    except:
        pass
    
record = BiopythonTranslator().translate_record(region_of_interest)
record = record.crop((b_start_p, b_end_p))
final_record["isescan"] = record

for feature in record.features:
    if feature.label == "ISEScan":
        feature.label = None

record.plot(ax = axs[1], x_lim=(b_start_p, b_end_p))
axs[1].annotate("ISEScan", (b_start_p - padding, 0), rotation = 90, annotation_clip=False)

# ax 2
window_size = 50
xx = list(range(b_start_p, b_end_p))
yy = gen_local_gc_content(region_of_interest.seq, window_size)
yy = [0] * int(window_size/2) + yy + [0] * int(window_size/2)
yy = yy[b_start_p: b_end_p]

axs[2].fill_between(xx, yy, alpha=0.3)
axs[2].set_ylim(bottom = 0)
axs[2].set_ylabel("GC(%)")

for i in axs[2].spines:
    axs[2].spines[i].set_visible(False)

for ax in axs:
    ax.set_xlim(b_start_p, b_end_p)

## Final
fig.suptitle(title)
plt.tight_layout()
plt.show()