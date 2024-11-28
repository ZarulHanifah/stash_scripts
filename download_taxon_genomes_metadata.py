#!/usr/bin/env python    

import sys
import os
import pathlib
import argparse    
import subprocess    
import json
import zipfile
import shutil
import pandas as pd
from tqdm import tqdm

# Function to run curl command and fetch taxon metadata    
def fetch_taxon_genome_details(taxon, sp_reps_only_bool=True):    
    sp_reps_only_bool2 = "true" if sp_reps_only_bool else "false"
    url = f'https://gtdb-api.ecogenomic.org/taxon/{taxon}/genomes-detail?sp_reps_only={sp_reps_only_bool2}'    
    result = subprocess.run(['curl', '-X', 'GET', url, '-H', 'accept: application/json'], capture_output=True, text=True)
    result = json.loads(result.stdout)
    result = result["rows"]
    return result

def fetch_genome_metadata(accession):
    url = f'https://gtdb-api.ecogenomic.org/genome/{accession}/card'
    result = subprocess.run(['curl', '-X', 'GET', url, '-H', 'accept: application/json'], capture_output=True, text=True)
    result = json.loads(result.stdout)
    return result

def parse_genome_metadata(gm):
    
    trna_aa_count    = gm["metadata_nucleotide"]["trna_aa_count"]
    contig_count     = gm["metadata_nucleotide"]["contig_count"]
    n50_contigs      = gm["metadata_nucleotide"]["n50_contigs"]
    longest_contig   = gm["metadata_nucleotide"]["longest_contig"]
    scaffold_count   = gm["metadata_nucleotide"]["scaffold_count"]
    n50_scaffolds    = gm["metadata_nucleotide"]["n50_scaffolds"]
    longest_scaffold = gm["metadata_nucleotide"]["longest_scaffold"]
    genome_size      = gm["metadata_nucleotide"]["genome_size"]
    gc_percentage    = gm["metadata_nucleotide"]["gc_percentage"]
    ambiguous_bases  = gm["metadata_nucleotide"]["ambiguous_bases"]

    checkm_completeness         = gm["metadata_gene"]["checkm_completeness"]
    checkm_contamination        = gm["metadata_gene"]["checkm_contamination"]
    checkm_strain_heterogeneity = gm["metadata_gene"]["checkm_strain_heterogeneity"]
    checkm2_completeness        = gm["metadata_gene"]["checkm2_completeness"]
    checkm2_contamination       = gm["metadata_gene"]["checkm2_contamination"]
    checkm2_model               = gm["metadata_gene"]["checkm2_model"]
    lsu_5s_count                = gm["metadata_gene"]["lsu_5s_count"]
    ssu_count                   = gm["metadata_gene"]["ssu_count"]
    lsu_23s_count               = gm["metadata_gene"]["lsu_23s_count"]
    protein_count               = gm["metadata_gene"]["protein_count"]
    coding_density              = gm["metadata_gene"]["coding_density"]
    
    ncbi_genbank_assembly_accession = gm["metadata_ncbi"]["ncbi_genbank_assembly_accession"]
    ncbi_strain_identifiers         = gm["metadata_ncbi"]["ncbi_strain_identifiers"]
    ncbi_assembly_level             = gm["metadata_ncbi"]["ncbi_assembly_level"]
    ncbi_assembly_name              = gm["metadata_ncbi"]["ncbi_assembly_name"]
    ncbi_assembly_type              = gm["metadata_ncbi"]["ncbi_assembly_type"]
    ncbi_bioproject                 = gm["metadata_ncbi"]["ncbi_bioproject"]
    ncbi_biosample                  = gm["metadata_ncbi"]["ncbi_biosample"]
    ncbi_country                    = gm["metadata_ncbi"]["ncbi_country"]
    ncbi_date                       = gm["metadata_ncbi"]["ncbi_date"]
    ncbi_genome_category            = gm["metadata_ncbi"]["ncbi_genome_category"]
    ncbi_genome_representation      = gm["metadata_ncbi"]["ncbi_genome_representation"]
    ncbi_isolate                    = gm["metadata_ncbi"]["ncbi_isolate"]
    ncbi_isolation_source           = gm["metadata_ncbi"]["ncbi_isolation_source"]
    ncbi_lat_lon                    = gm["metadata_ncbi"]["ncbi_lat_lon"]
    ncbi_molecule_count             = gm["metadata_ncbi"]["ncbi_molecule_count"]
    ncbi_cds_count                  = gm["metadata_ncbi"]["ncbi_cds_count"]
    ncbi_refseq_category            = gm["metadata_ncbi"]["ncbi_refseq_category"]
    ncbi_seq_rel_date               = gm["metadata_ncbi"]["ncbi_seq_rel_date"]
    ncbi_spanned_gaps               = gm["metadata_ncbi"]["ncbi_spanned_gaps"]
    ncbi_species_taxid              = gm["metadata_ncbi"]["ncbi_species_taxid"]
    ncbi_ssu_count                  = gm["metadata_ncbi"]["ncbi_ssu_count"]
    ncbi_submitter                  = gm["metadata_ncbi"]["ncbi_submitter"]
    ncbi_taxid                      = gm["metadata_ncbi"]["ncbi_taxid"]
    ncbi_total_gap_length           = gm["metadata_ncbi"]["ncbi_total_gap_length"]
    ncbi_translation_table          = gm["metadata_ncbi"]["ncbi_translation_table"]
    ncbi_trna_count                 = gm["metadata_ncbi"]["ncbi_trna_count"]
    ncbi_unspanned_gaps             = gm["metadata_ncbi"]["ncbi_unspanned_gaps"]
    ncbi_version_status             = gm["metadata_ncbi"]["ncbi_version_status"]
    ncbi_wgs_master                 = gm["metadata_ncbi"]["ncbi_wgs_master"]

    gtdbTypeDesignation        = gm["metadata_type_material"]["gtdbTypeDesignation"]
    gtdbTypeDesignationSources = gm["metadata_type_material"]["gtdbTypeDesignationSources"]
    lpsnTypeDesignation        = gm["metadata_type_material"]["lpsnTypeDesignation"]
    dsmzTypeDesignation        = gm["metadata_type_material"]["dsmzTypeDesignation"]
    lpsnPriorityYear           = gm["metadata_type_material"]["lpsnPriorityYear"]
    gtdbTypeSpeciesOfGenus     = gm["metadata_type_material"]["gtdbTypeSpeciesOfGenus"]

    ncbi_taxonomy                  = gm["metadataTaxonomy"]["ncbi_taxonomy"]
    ncbi_taxonomy_unfiltered       = gm["metadataTaxonomy"]["ncbi_taxonomy_unfiltered"]
    gtdb_representative            = gm["metadataTaxonomy"]["gtdb_representative"]
    gtdb_genome_representative     = gm["metadataTaxonomy"]["gtdb_genome_representative"]
    ncbi_type_material_designation = gm["metadataTaxonomy"]["ncbi_type_material_designation"]
    gtdbDomain                     = gm["metadataTaxonomy"]["gtdbDomain"]
    gtdbPhylum                     = gm["metadataTaxonomy"]["gtdbPhylum"]
    gtdbClass                      = gm["metadataTaxonomy"]["gtdbClass"]
    gtdbOrder                      = gm["metadataTaxonomy"]["gtdbOrder"]
    gtdbFamily                     = gm["metadataTaxonomy"]["gtdbFamily"]
    gtdbGenus                      = gm["metadataTaxonomy"]["gtdbGenus"]
    gtdbSpecies                    = gm["metadataTaxonomy"]["gtdbSpecies"]
    
    result = [trna_aa_count, contig_count, n50_contigs, longest_contig, scaffold_count, n50_scaffolds, longest_scaffold, genome_size, gc_percentage, ambiguous_bases, checkm_completeness, checkm_contamination, checkm_strain_heterogeneity, checkm2_completeness, checkm2_contamination, checkm2_model, lsu_5s_count, ssu_count, lsu_23s_count, protein_count, coding_density, ncbi_genbank_assembly_accession, ncbi_strain_identifiers, ncbi_assembly_level, ncbi_assembly_name, ncbi_assembly_type, ncbi_bioproject, ncbi_biosample, ncbi_country, ncbi_date, ncbi_genome_category, ncbi_genome_representation, ncbi_isolate, ncbi_isolation_source, ncbi_lat_lon, ncbi_molecule_count, ncbi_cds_count, ncbi_refseq_category, ncbi_seq_rel_date, ncbi_spanned_gaps, ncbi_species_taxid, ncbi_ssu_count, ncbi_submitter, ncbi_taxid, ncbi_total_gap_length, ncbi_translation_table, ncbi_trna_count, ncbi_unspanned_gaps, ncbi_version_status, ncbi_wgs_master, gtdbTypeDesignation, gtdbTypeDesignationSources, lpsnTypeDesignation, dsmzTypeDesignation, lpsnPriorityYear, gtdbTypeSpeciesOfGenus, ncbi_taxonomy, ncbi_taxonomy_unfiltered, gtdb_representative, gtdb_genome_representative, ncbi_type_material_designation, gtdbDomain, gtdbPhylum, gtdbClass, gtdbOrder, gtdbFamily, gtdbGenus, gtdbSpecies]
    
    return result

# Parse command line arguments    
def parse_arguments():    
    parser = argparse.ArgumentParser(description="Fetch metadata for a given taxon")    
    parser.add_argument('-s', '--taxon', required=True, help="Taxon string to search (e.g., 'c__Bog-38')")    
    parser.add_argument('--not-sp-reps-only', action='store_true', default=False,
            help="If specified, include non-species representatives (default is species representatives only)")
    parser.add_argument('-o', '--output', type=str, help="Output file to save the the results (optional, default is stdout)")
    parser.add_argument('--download-genomes', action='store_true', default=False,
            help="If specified, download genomes for each accession")
    parser.add_argument('--genome-dir', type=str, help="Directory to save downloaded genomes (required if --download-genomes is specified)")

    return parser.parse_args()    

# Main function    
def main():    
    # Parse arguments     
    args = parse_arguments()
    if args.download_genomes and not args.genome_dir:
        raise ValueError("You must specify --genome-dir when --download-genomes is used.")
    if args.download_genomes:
        args.output = f"{args.genome_dir}/metadata.tsv"
        os.makedirs(args.genome_dir, exist_ok=True)

    taxon = args.taxon    
    sp_reps_only_bool = not args.not_sp_reps_only  # Default is True, unless --not-sp-reps-only is specified
    
    # Fetch metadata for the given taxon    
    metadata = fetch_taxon_genome_details(taxon, sp_reps_only_bool)    
    # print(metadata["rows"])
    
    # Print the metadata    
    # print(json.dumps(metadata, indent=4))    
    
    columns_names = ["trna_aa_count", "contig_count", "n50_contigs", "longest_contig", "scaffold_count", "n50_scaffolds", "longest_scaffold", "genome_size", "gc_percentage", "ambiguous_bases", "checkm_completeness", "checkm_contamination", "checkm_strain_heterogeneity", "checkm2_completeness", "checkm2_contamination", "checkm2_model", "lsu_5s_count", "ssu_count", "lsu_23s_count", "protein_count", "coding_density", "ncbi_genbank_assembly_accession", "ncbi_strain_identifiers", "ncbi_assembly_level", "ncbi_assembly_name", "ncbi_assembly_type", "ncbi_bioproject", "ncbi_biosample", "ncbi_country", "ncbi_date", "ncbi_genome_category", "ncbi_genome_representation", "ncbi_isolate", "ncbi_isolation_source", "ncbi_lat_lon", "ncbi_molecule_count", "ncbi_cds_count", "ncbi_refseq_category", "ncbi_seq_rel_date", "ncbi_spanned_gaps", "ncbi_species_taxid", "ncbi_ssu_count", "ncbi_submitter", "ncbi_taxid", "ncbi_total_gap_length", "ncbi_translation_table", "ncbi_trna_count", "ncbi_unspanned_gaps", "ncbi_version_status", "ncbi_wgs_master", "gtdbTypeDesignation", "gtdbTypeDesignationSources", "lpsnTypeDesignation", "dsmzTypeDesignation", "lpsnPriorityYear", "gtdbTypeSpeciesOfGenus", "ncbi_taxonomy", "ncbi_taxonomy_unfiltered", "gtdb_representative", "gtdb_genome_representative", "ncbi_type_material_designation", "gtdbDomain", "gtdbPhylum", "gtdbClass", "gtdbOrder", "gtdbFamily", "gtdbGenus", "gtdbSpecies"]
    index_names = []
    df = []

    for genome in metadata:
        acc = genome["gid"]
        index_names.append(acc)
        gm = fetch_genome_metadata(acc)
        gm = parse_genome_metadata(gm)
        df.append(gm)

    if args.download_genomes:
        for genome in tqdm(metadata, desc="Downloading genomes"):
            acc = genome["gid"]
            genome_path = pathlib.Path(args.genome_dir) / f"{acc}.fasta"
            # genome_path = os.path.join(args.genome_dir, f"{acc}.fasta")
            if genome_path.exists():
                print(f"Skipping download for {acc}, {genome_path} already exists.", file=sys.stderr)
                continue
            subprocess.run(['rm', '-rf', 'ncbi_dataset.zip'], check=True)
            subprocess.run(['datasets', 'download', 'genome', 'accession', acc,
                            '--include', 'genome', '--no-progressbar', '--filename', 'ncbi_dataset.zip'])

            with zipfile.ZipFile('ncbi_dataset.zip', 'r') as zip_ref:
                fna_files = [f for f in zip_ref.namelist() if f.endswith(".fna")]
                if not fna_files:
                    raise ValueError(f"No .fna file found in ncbi_dataset.zip for accession {acc}.")
                # Extract only the .fna file and write to the target path
                with zip_ref.open(fna_files[0]) as source, open(genome_path, 'wb') as target:
                    target.write(source.read())

            subprocess.run(['rm', '-rf', 'ncbi_dataset.zip'], check=True)
            # print(f"Genome {acc} saved to {genome_path}.")

    df = pd.DataFrame(df, index=index_names, columns=columns_names)

    if args.output:
        df.to_csv(args.output, sep="\t")
    else:
        df.to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":    
    main()

