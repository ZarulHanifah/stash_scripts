#!/usr/bin/env python3

from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt
from BCBio import GFF
import pandas as pd
import numpy as np

class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)() # retain local pointer to value
        return value                     # faster to return than dict lookup

dram_df = pd.read_csv("annotations.tsv", sep = "\t", index_col = 0)
dram_df.sort_index(inplace = True)

target_ko = "K16157"
target_list = []

for g in dram_df.index.tolist():
    if dram_df.loc[g, "ko_id"] == target_ko:
        target_list.append(g)
        

gb_path = "/home/zarul/Zarul/Peatlands/status_nanopore/Methylovirgula/bakta/all_annotations/GCF_019343105.1/GCF_019343105.1.gff3"
gff_records = []

with open(gb_path) as in_handle:
    for rec in GFF.parse(in_handle):
        gff_records.append(rec)

left_shift = 3
right_shift = 7

hits = Vividict()

def get_slice_given_id(feature_list, gene_id):
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
            return f


# gen list of gene range
for r in gff_records:
    contig_id = r.id
    
    len_features = len(r.features)
    list_features = [f.id for f in r.features]
    
    for f in r.features:
        if f.id in target_list:
            idx = target_list.index(f.id)
            
            hits[idx]["contig"] = contig_id
            hits[idx]["gene_ids"] = get_slice_given_id(list_features, f.id)
            ids_range = get_slice_given_id(list_features, f.id)
            # hits[idx]["genes"] = [get_record(gff_records, contig_id, i) for i in ids_range]
            hits[idx]["start"] = int(get_record(gff_records, contig_id, ids_range[0]).location.start)
            hits[idx]["end"] = int(get_record(gff_records, contig_id, ids_range[-1]).location.end)
    

fig, ax = plt.subplots(figsize = (15, 6), dpi = 300)

for h in list(hits.keys()):
    region_of_interest = [i for i in gff_records if i.id == hits[h]["contig"]][0]
    seq = region_of_interest._seq
    region_of_interest.features = [f for f in region_of_interest.features if f.id in hits[h]["gene_ids"]]
    
    for i in region_of_interest.features:
        try:
            gene_id = i.qualifiers["ID"][0]
            product = dram_df.loc[gene_id, "ko_id"]
            if product is not np.nan:
                i.qualifiers["product"] = f"{gene_id}:\n{product}"
        except:
            pass
        
    record = BiopythonTranslator().translate_record(region_of_interest)
    record = record.crop((hits[h]["start"], hits[h]["end"]))
    
    # for feature in record.features:
    #     feature.label = None
        
    record.plot(ax = ax)