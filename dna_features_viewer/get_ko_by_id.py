#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 14:44:37 2023

@author: ahbui
"""

from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt
from BCBio import GFF
import pandas as pd
import numpy as np
import re

pd.set_option("display.max_colwidth", None)

class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)() # retain local pointer to value
        return value                     # faster to return than dict lookup

# tg = "GCF_019343105.1"     # Gwak, HY1
# tg = "GCF_004135935.1"     # ligni
tg = "metabat.73"

left_shift = 5
right_shift = 5
orient = 1

# target_ko = "K16157"     # mmoX
# target_ko = "K14063"     # AraC transcriptional regulator
target_ko = "PF09407.13"   # AbiEi anti-toxin
target_ko = "K01601"       # rubisco large chain

dram_df = pd.read_csv(f"../dram_faa/{tg}/annotations.tsv", sep = "\t", index_col = 0)
dram_df.sort_index(inplace = True)

gb_path = f"../bakta/{tg}/{tg}.gff3"

naming = pd.read_csv("naming_table.tsv", sep = "\t", index_col = 0).to_dict(orient = "index")

target_list = []

def get_slice_given_id(feature_list, gene_id, left_shift, right_shift):
    f_index = feature_list.index(gene_id)
    start = max(f_index - left_shift, 0)
    end = min(f_index + right_shift + 1, len_features)
    ids_range = list_features[start:end]
    return ids_range

# given a gene id, get the record
def get_record(gff_records, contig_id, gene_id):
    list_of_contigs = [gff_records[i].id for i in range(len(gff_records))]
    contig_idx = list_of_contigs.index(contig_id)
    
    t = gff_records[contig_idx]
    for f in t.features:
        if f.id == gene_id:
            if f.type == "region":
                return
            return f

def clean_pfam_product(p):
    product = re.sub(".*\[", "", p)
    product = re.sub("\].*", "", product)
    product = re.sub("\].*", "", product)
    return product

def last_resort_get_product(gene_id):
    with open(gb_path, "r") as f:
        ls = f.readlines()
    for l in ls:
        if gene_id in l:
            final = re.sub(".*Name=", "", l)
            final = re.sub(";.*", "", final)
            final = final.strip()
            return f"bakta: {final}"
        
class MyCustomTranslator(BiopythonTranslator):
    def __init__(self, target_id, features_filters=(), features_properties=None):
        self.target_id = target_id
        self.features_filters = features_filters
        self.features_properties = features_properties
        
    def compute_feature_color(self, feature):
        if feature.id == self.target_id:
            return "green"
        # else:
        #     return "purple"
        try:
            rank = dram_df.loc[feature.id, "rank"]
            if rank == "C":
                product = dram_df.loc[feature.id, "ko_id"]
            elif rank == "D":
                product = dram_df.loc[feature.id, "pfam_hits"]
                product = clean_pfam_product(product)
            color = naming[product]["color"]
            return color
        except:
            return "#0000ff"

for gene_id in dram_df.index.tolist():
    
    if dram_df.loc[gene_id, "ko_id"] == target_ko:
        target_list.append(gene_id)
    elif isinstance(dram_df.loc[gene_id, "pfam_hits"], str):
        pfam_id = clean_pfam_product(dram_df.loc[gene_id, "pfam_hits"])
        
        if pfam_id == target_ko:
            target_list.append(gene_id)

gff_records = []
with open(gb_path) as in_handle:
    for rec in GFF.parse(in_handle):
        gff_records.append(rec)

hits = Vividict()

# gen list of gene range
for r in gff_records:
    contig_id = r.id
    
    len_features = len(r.features)
    list_features = [f.id for f in r.features if f.type not in ["region", "gap"]]
    
    for f in r.features:
        if f.id in target_list:
            idx = target_list.index(f.id)
            
            hits[idx]["contig"] = contig_id
            print(contig_id)
            hits[idx]["anchor"] = f.id
            hits[idx]["strand"] = f.strand
            ids_range = get_slice_given_id(list_features, f.id, left_shift, right_shift)
            # if f.strand == orient:
            # else:
                # ids_range = get_slice_given_id(list_features, f.id, right_shift, left_shift)
            
            hits[idx]["gene_ids"] = ids_range
            # hits[idx]["genes"] = [get_record(gff_records, contig_id, i) for i in ids_range]
            hits[idx]["start"] = int(get_record(gff_records, contig_id, ids_range[0]).location.start)
            hits[idx]["end"] = int(get_record(gff_records, contig_id, ids_range[-1]).location.end)
    
nrows = len(list(hits.keys()))
fig, axs = plt.subplots(nrows = nrows, figsize = (15, 5 * nrows), dpi = 100)

df = []
for ax_idx, h in enumerate(list(hits.keys())):
    region_of_interest = [i for i in gff_records if i.id == hits[h]["contig"]][0]
    seq = region_of_interest._seq
    region_of_interest.features = [f for f in region_of_interest.features if f.id in hits[h]["gene_ids"]]
    
    for i in region_of_interest.features:
        rank, product, strand, short, long = np.nan, np.nan, np.nan, np.nan, np.nan
        
        gene_id = i.qualifiers["ID"][0]
        try:
            rank = dram_df.loc[gene_id, "rank"]
        except:
            pass
        
        if rank == "C":
            product = dram_df.loc[gene_id, "ko_id"]
            long = dram_df.loc[gene_id, "kegg_hit"]
        elif rank == "D":
            product = dram_df.loc[gene_id, "pfam_hits"]
            product = clean_pfam_product(product)
            long = dram_df.loc[gene_id, "pfam_hits"]
        elif rank == "E":
            product = "hypothetical"
        else:
            product = last_resort_get_product(gene_id)

        try:
            short = naming[product]["short"]
            # long= naming[product]["long"]
            if short is not np.nan:
                i.qualifiers["gene"] = f"{short}"
        except:
            i.qualifiers["gene"] = re.sub(".*_", "_", gene_id)
        
        df.append([gene_id, rank, product, short, long])
    
    record = MyCustomTranslator(hits[h]["anchor"]).translate_record(region_of_interest)
    record = record.crop((hits[h]["start"], hits[h]["end"]))
    
    # for feature in record.features:
    #     feature.label = None
    if nrows == 1:
        ax = axs
    else:
        ax = axs[ax_idx]
    record.plot(ax = ax, annotate_inline = False)
    
    # reorient if you want
    if hits[h]["strand"] != orient:
        x_start, x_end = ax.get_xlim()
        ax.set_xlim(x_end, x_start)

if orient == -1:
    df = list(reversed(df))


df = pd.DataFrame(df)
df.columns = ["Gene ID", "DRAM rank", "product", "short", "long"]
df.set_index("Gene ID", inplace = True)
df.drop_duplicates(inplace = True)
print(df)

df.to_csv("summary.csv", sep = "\t")