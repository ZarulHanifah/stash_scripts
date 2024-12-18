import os
import re
import subprocess
import pandas as pd
import dask.dataframe as dd
from Bio.KEGG import REST
import ipywidgets as widgets
from IPython.display import display
from multiprocessing import Process, Manager
from concurrent.futures import ProcessPoolExecutor, as_completed
import plotly.express as px
from plotly.subplots import make_subplots
from IPython.display import clear_output

import dask.dataframe as dd
from tqdm.notebook import tqdm

class KEGGModuleDashboard:
    def __init__(self):
        # Step 1: Fetch KEGG modules and prepare data
        self.module_df = self.fetch_kegg_modules()
        self.modules_formatted_names = (
            self.module_df["module"] + ": " + self.module_df["name"]
        ).tolist()

        # Step 2: Define widgets
        self.filter_text = widgets.Text(
            value="",
            description="Filter modules",
            placeholder="Type to filter..."
        )

        self.modules_select_w = widgets.Select(
            options=self.modules_formatted_names,
            description="Modules list",
            rows=20,
            layout=widgets.Layout(width="450px")
        )

        self.module_info_output = widgets.Output()

        # Step 3: Attach widget behaviors
        self.filter_text.observe(self.filter_modules, names="value")
        self.modules_select_w.observe(self.display_module_info, names="value")

        # Step 4: Create the dashboard layout
        self.dashboard = widgets.HBox(
            [widgets.VBox([self.filter_text, self.modules_select_w]), self.module_info_output]
        )

    def fetch_kegg_modules(self):
        """Fetch the list of KEGG modules and return as a DataFrame."""
        module_list_data = REST.kegg_list("module").read()
        return pd.DataFrame(
            [line.split("\t") for line in module_list_data.strip().split("\n")],
            columns=["module", "name"]
        )

    def filter_modules(self, change):
        """Update the options in the module select widget based on filter text."""
        search_text = change["new"]
        if search_text == "":
            filtered_modules = self.modules_formatted_names
        else:
            filtered_modules = [
                name for name in self.modules_formatted_names
                if search_text.lower() in name.lower()
            ]
        self.modules_select_w.options = filtered_modules

    def display_module_info(self, change):
        """Fetch and display module information for the selected module."""
        self.module_info_output.clear_output()
        selected_module = change["new"].split(":")[0]  # Extract module ID
        if selected_module:
            with self.module_info_output:
                module_info = REST.kegg_get(f"module:{selected_module}").readlines()
                for l in module_info:
                    if "NAME" in l:
                        name = l
                        break
                name = re.sub("NAME *", "", name)
                print(f"{selected_module}: {name}")
                module_info_description = "".join(module_info)
                print(module_info_description)

    def display(self):
        """Display the dashboard."""
        display(self.dashboard)

class KOHeatmapDashboard:
    def __init__(self, dram_dir, ko_list_path, ko_synonyms_path, summary_all):
        # Load KO list and synonyms
        self.genomes_annotations = {}
        self.dram_dir = dram_dir
        self.ko_list_df = self.load_ko_list(ko_list_path, ko_synonyms_path)
        self.ko_list = (self.ko_list_df.index.compute() + ": " + self.ko_list_df["gene_name"].compute()).to_dict()
        self.ko_formatted_names = (
            (self.ko_list_df.index.compute() + ": " + self.ko_list_df["gene_name"].compute())
        ).tolist()
        self.summary_all = summary_all
        
        # Custom color scheme
        self.colors = ["#e8e8e8", "#df0024", "#0085c7"]

        self.read_csv_dtype = {'fasta': "object",
            'rank': "object",
            'ko_id': "object",
            'kegg_hit': "object",
            'peptidase_id': "object",
            'peptidase_family': "object",
            'peptidase_hit': "object",
            'peptidase_RBH': "object",
            'peptidase_identity': "float64",
            'peptidase_bitScore': "float64",
            'peptidase_eVal': "float64",
            'pfam_hits': "object",
            'cazy_ids': "object",
            'cazy_hits': "object",
            'cazy_subfam_ec': "object",
            'cazy_best_hit': "object",
            'vogdb_id': "object",
            'vogdb_hits': "object",
            'vogdb_categories': "object",
            'heme_regulatory_motif_count': "float64"}

        # Global matrix for persistence
        self.df_count_ko = pd.DataFrame()

        # Widgets
        self.filter_text = widgets.Text(
            value="",
            description="Filter KO list",
            placeholder="Type to filter..."
        )
        self.ko_formatted_names_w = widgets.Select(
            options=self.ko_formatted_names,
            description="KO list",
            rows=40,
            layout=widgets.Layout(width='450px')
        )
        self.genomes_input = widgets.TagsInput(
            value=[],
            description="Genomes",
            allow_duplicates=False
        )
        self.ko_input = widgets.TagsInput(
            value=[],
            description="KEGG KOs",
            allow_duplicates=False
        )
        self.heatmap_output = widgets.Output()
        self.reset_button = widgets.Button(description="Reset Heatmap")
        
        # Attach behavior
        self.genomes_input.observe(self.update_heatmap, names="value")
        self.ko_input.observe(self.update_heatmap, names="value")
        self.filter_text.observe(self.filter_ko_list, names="value")
        self.ko_formatted_names_w.observe(self.add_to_ko_input, names="value")
        self.ko_formatted_names_w.value = None
        self.reset_button.on_click(self.reset_table_and_replot)
        self.updating_heatmap = False  # Add a flag to prevent duplicate updates

    def load_ko_list(self, ko_list_path, ko_synonyms_path):
        """Load and format the KO list with synonyms."""
        ko_list_df = dd.read_csv(ko_list_path, sep="\t")
        ko_list_df = ko_list_df.set_index(ko_list_df.columns[0])

        ko_list_df = ko_list_df.loc[:, ["definition"]]
        ko_list_df.columns = ["gene_name"]


        ko_synonyms = dd.read_csv(ko_synonyms_path, sep="\t", header=None) # , names=["synonyms"]
        ko_synonyms = ko_synonyms.set_index(ko_synonyms.columns[0])
        ko_synonyms = ko_synonyms.fillna("")
        ko_synonyms.columns = ["synonyms"]

        ko_list_df = ko_list_df.merge(ko_synonyms, left_index=True, right_index=True)
        ko_list_df["gene_name"] = ko_list_df["synonyms"] + ", " + ko_list_df["gene_name"]
        ko_list_df = ko_list_df.drop(columns=["synonyms"])
        ko_list_df.index.name = "KO"
        return ko_list_df

    def read_df(self, sample):
        """Optimized function to read sample data."""
        if sample in self.genomes_annotations:
            return self.genomes_annotations[sample]
        annotation_path = os.path.join(self.dram_dir, f"{sample}/annotations.tsv")
        annot_df = pd.read_csv(annotation_path, sep="\t", dtype=self.read_csv_dtype)
        self.genomes_annotations[sample] = annot_df  # Store in memory
        return annot_df

    def count_sample_ko(self, annot_df, ko):
        """Optimized function to count KO occurrences."""
        return annot_df["ko_id"].eq(ko).sum()

    def process_genomes_input(self, genomes_input_value):
        processed_genomes = []
        for genome in genomes_input_value:
            genome = re.sub(r"[\'\[\],]", " ", genome)
            if " " in genome:  # Check for whitespace
                new_genomes = genome.split()
                processed_genomes.extend(new_genomes)  # Split and add the items in order
            else:
                processed_genomes.append(genome)
        return processed_genomes

    def process_ko_input(self, ko_input_value):
        processed_kos = []
        for ko in self.ko_input.value:
            ko = re.sub(r"[\+\-,()]", " ", ko)
            if " " in ko:
                kos = ko.split()
                kos = [re.sub(r"['\[\]\+\-,()]", " ", ko) for ko in kos]
                kos = [ko.strip() for ko in kos]
                kos = [k for k in kos if re.match("^K[0-9]*", k)]
                kos = remove_duplicates(kos)
                processed_kos.extend(kos)
            else:
                processed_kos.append(ko)
        return processed_kos

    def update_matrix(self, df_count_ko, genomes, genes):
        """Update the KO-genome occurrence matrix."""
        new_columns = [self.ko_list[gene] for gene in genes]
        df_count_ko = df_count_ko.reindex(columns=new_columns)

        for gene in genes:
            if self.ko_list[gene] not in df_count_ko.columns:
                df_count_ko[self.ko_list[gene]] = 0

        for genome in genomes:
            if genome not in df_count_ko.index:
                df_count_ko.loc[genome] = 0

        for genome in genomes:
            for gene in genes:
                gene_name = self.ko_list[gene]
                counts = self.count_sample_ko(self.read_df(genome), gene)
                df_count_ko.at[genome, gene_name] = counts

        df_count_ko = df_count_ko.loc[genomes, :]
        return df_count_ko

    def add_hover_info(self, fig, summary_all, df_count_ko, selected_genomes):
        # Extract classification and additional metadata
        genome_metadata = {
            genome: {
                "Classification": summary_all.loc[genome, "classification"] if genome in summary_all.index else "N/A",
                "Completeness": summary_all.loc[genome, "Completeness"] if genome in summary_all.index else "N/A",
                "Contamination": summary_all.loc[genome, "Contamination"] if genome in summary_all.index else "N/A",
                "number_of_contigs": summary_all.loc[genome, "number_of_contigs"] if genome in summary_all.index else "N/A",
            }
            for genome in selected_genomes
        }

        # Generate customdata for the heatmap
        hover_metadata = [
            [
                [
                    genome_metadata.get(row, {}).get("Classification", "N/A"),
                    genome_metadata.get(row, {}).get("Completeness", "N/A"),
                    genome_metadata.get(row, {}).get("Contamination", "N/A"),
                    genome_metadata.get(row, {}).get("number_of_contigs", "N/A")
                ]
                for _ in range(df_count_ko.shape[1])
            ]
            for row in df_count_ko.index
        ]

        fig.update_traces(
            customdata=hover_metadata,
            hovertemplate=(
                "Genome: %{y}<br>"
                "KO: %{x}<br>"
                "Count: %{z}<br>"
                "Classification: %{customdata[0]}<br>"
                "Completeness: %{customdata[1]}<br>"
                "Contamination: %{customdata[2]}<br>"
                "Number of Contigs: %{customdata[3]}"
            )
        )

    def update_heatmap(self, change=None):
        """Update the heatmap visualization."""
        if self.updating_heatmap:
            return  # Skip if already updating
        
        self.updating_heatmap = True  # Set the flag
        with self.heatmap_output:
            clear_output(wait=True)
        
            self.genomes_input.value = self.process_genomes_input(self.genomes_input.value)
            selected_genomes = self.genomes_input.value
            
            # Read annotations for each genome
            for genome in selected_genomes:
                if genome not in self.genomes_annotations:
                    self.read_df(genome)

            self.ko_input.value = self.process_ko_input(self.ko_input.value)
            selected_kos = self.ko_input.value
        
            self.df_count_ko = self.update_matrix(self.df_count_ko, selected_genomes, selected_kos)
        
            num_rows, num_cols = self.df_count_ko.shape
            cell_size = 30  # Size per cell in pixels
            min_width = 400  # Minimum width for heatmap
            min_height = 300  # Minimum height for heatmap
        
            # Dynamically adjust the figure dimensions
            heatmap_width = max(num_cols * cell_size, min_width)
            heatmap_height = max(num_rows * cell_size, min_height)
            heatmap_height += 300

            fig_dims = {  #   width, height, col
                "completeness_fig": [60, heatmap_height, 1],
                "contamination_fig": [60, heatmap_height, 2],
                "coverage_fig": [60, heatmap_height, 3],
                "n_contigs_fig": [60, heatmap_height, 4],
                "heatmap_fig": [heatmap_width, heatmap_height, 5]
            }

            margin_dict = dict(l=0, t=50, r=0, b=0)

            completeness_fig = px.bar(self.summary_all.loc[selected_genomes, :],
                y=self.summary_all.loc[selected_genomes, :].index, x="Completeness", orientation="h")

            contamination_fig = px.bar(self.summary_all.loc[selected_genomes, :],
                y=self.summary_all.loc[selected_genomes, :].index, x="Contamination", orientation="h")

            coverage_fig = px.bar(self.summary_all.loc[selected_genomes, :],
                y=self.summary_all.loc[selected_genomes, :].index, x="average_coverage", orientation="h")

            n_contigs_fig = px.bar(self.summary_all.loc[selected_genomes, :],
                y=self.summary_all.loc[selected_genomes, :].index, x="number_of_contigs", orientation="h")

            # Create the heatmap plot
            heatmap_fig = px.imshow(self.df_count_ko, aspect="auto", text_auto=True)
            heatmap_fig.update_layout(width=heatmap_width, height=heatmap_height, margin=margin_dict)

            self.add_hover_info(heatmap_fig, self.summary_all, self.df_count_ko, selected_genomes)

            column_widths = [fig_dims[f][0] for f in fig_dims]
            combined_fig = make_subplots(
                rows=1, cols=len(fig_dims), shared_yaxes=True, column_widths=column_widths, horizontal_spacing=0.005
            )

            # COMPLETENESS
            width, height, col = fig_dims["completeness_fig"]
            combined_fig.add_trace(completeness_fig.data[0], row=1, col=col)
            combined_fig.update_xaxes(range=[50,100], title=dict(text="Comp.(%)", font=dict(size=12)), row=1, col=col)
            
            # CONTAMINATION
            width, height, col = fig_dims["contamination_fig"]
            combined_fig.add_trace(contamination_fig.data[0], row=1, col=col)
            combined_fig.update_xaxes(range=[0,10], title=dict(text="Cont.(%)", font=dict(size=12)), row=1, col=col)


            # COVERAGE
            width, height, col = fig_dims["coverage_fig"]
            combined_fig.add_trace(coverage_fig.data[0], row=1, col=col)
            combined_fig.update_xaxes(title=dict(text="Cov (x)", font=dict(size=12)), row=1, col=col)

            # N_CONTIGS
            width, height, col = fig_dims["n_contigs_fig"]
            combined_fig.add_trace(n_contigs_fig.data[0], row=1, col=col)
            combined_fig.update_xaxes(title=dict(text="n contigs", font=dict(size=12)), row=1, col=col)

            
            width, height, col = fig_dims["heatmap_fig"]
            combined_fig.add_trace(heatmap_fig.data[0], row=1, col=col)

            tickvals = list(range(self.df_count_ko.shape[1]))
            xticktext = [re.sub(",.*", "", c) for c in self.df_count_ko.columns.tolist()]
            combined_fig.update_xaxes(tickangle=90, ticktext=xticktext, tickvals=tickvals, row=1, col=col)


            # Apply the color scale to the heatmap trace using coloraxis
            combined_fig.update_traces(coloraxis="coloraxis2", row=1, col=col)

            # Update coloraxes for the second subplot (heatmap)
            total_width = sum([fig_dims[f][0] for f in fig_dims])

            combined_fig.update_layout(
                coloraxis2=dict(
                    colorscale=self.colors, cmin=0, cmax=2, 
                    colorbar=dict(tickvals=[0, 1, 2], ticktext=["0", "1", "2+"]),
                ),
                coloraxis2_showscale=False,
                width=total_width, height=heatmap_height, font=dict(size=15), margin=margin_dict
            )
            combined_fig.update_yaxes(autorange="reversed")
            
            # combined_fig.show()
            combined_fig.show(config={"editable": True})
        self.updating_heatmap = False  # Reset the flag


    def reset_table_and_replot(self, change=None):
        """Reset the matrix and refresh the heatmap."""
        self.df_count_ko = pd.DataFrame()
        selected_genomes = self.genomes_input.value
        selected_kos = self.ko_input.value
        new_columns = [self.ko_list[gene] for gene in selected_kos]
        self.df_count_ko = self.df_count_ko.reindex(columns=new_columns)
        self.df_count_ko = self.update_matrix(self.df_count_ko, selected_genomes, selected_kos)
        self.update_heatmap(None)

    def filter_ko_list(self, change):
        """Filter KO list based on user input."""
        self.ko_formatted_names_w.value = None
        filter_value = self.filter_text.value.lower()
        filtered_options = [
            ko for ko in self.ko_formatted_names if filter_value in ko.lower()
        ]
        self.ko_formatted_names_w.options = filtered_options

        if not filter_value:
            self.ko_formatted_names_w.value = None

    def add_to_ko_input(self, change):
        """Add selected KO to the input list."""
        selected_item = self.ko_formatted_names_w.value
        if selected_item:
            selected_ko = selected_item.split(":")[0]
            if selected_ko not in self.ko_input.value:
                self.ko_input.value = self.ko_input.value + [selected_ko]

    def display(self):
        """Display the dashboard."""
        layout = widgets.VBox([
            widgets.HBox([
                self.heatmap_output,
                widgets.VBox(),
                widgets.VBox([self.filter_text, self.ko_formatted_names_w, self.reset_button])
            ]),
            widgets.HTML(value="<b>Selected Genomes </b>"),
            self.genomes_input,
            widgets.HTML(value="<b>Selected KEGG KO</b>"),
            self.ko_input
        ])
        display(layout)


# Usage
# dashboard = KOHeatmapDashboard()
# dashboard.display()

def remove_duplicates(input_list):
    seen = set()
    return [x for x in input_list if x not in seen and not seen.add(x)]

def does_my_genome_have_this_kegg_gene_with_grep2(dram_dir, genome, kegg_gene):
    # Construct the file path
    file_path = os.path.join(dram_dir, f"{genome}/annotations.tsv")
    
    # Use subprocess to run grep and capture output, looking for exact match in the 4th column
    result = subprocess.run(
        ["grep", "-P", "^([^\t]*\t){3}" + kegg_gene, file_path],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    
    # If grep found the gene, return True, otherwise False
    return result.returncode == 0  # grep returns 0 if found, non-zero if not

def filter_genomes_for_kegg_gene(dram_dir, genomes, gene, max_workers=6):
    """Filter genomes to check if they contain a specific KEGG gene using multiprocessing.
    
    Args:
        genomes (list): List of genomes
        gene (str): KEGG ID of gene
        max_workers (int): Number of workers
    """

    filtered_genomes = []
    
    # Use ProcessPoolExecutor for multiprocessing
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit tasks to the pool
        futures = {
            executor.submit(does_my_genome_have_this_kegg_gene_with_grep2, dram_dir, g, gene): g
            for g in genomes
        }

        # Use tqdm to track progress
        for future in tqdm(as_completed(futures), total=len(futures), desc=f"Filtering genomes for {gene}"):
            genome = futures[future]
            try:
                # Check the result of each task
                if future.result():
                    filtered_genomes.append(genome)
            except Exception as e:
                print(f"Error processing genome {genome}: {e}")

    return filtered_genomes

class quickSearchKObyGeneName:
    def __init__(self, ko_list_path, ko_synonyms_path, show_genes=False):
        """Initialize the quickSearchKObyGeneName class.

        Args:
            ko_list_path (str): Path to the KO list file.
            ko_synonyms_path (str): Path to the KO synonyms file.
            show_genes (bool): Show genes info or not?

        Raises:
            FileNotFoundError: If the specified file paths do not exist.
        """
        # Validate file paths
        self.validate_file_path(ko_list_path, "KO list file")
        self.validate_file_path(ko_synonyms_path, "KO synonyms file")

        # Load KO list and synonyms
        self.ko_list_df = self.load_ko_list(ko_list_path, ko_synonyms_path)
        self.ko_list = (self.ko_list_df.index + ": " + self.ko_list_df["gene_name"]).to_dict()
        self.ko_formatted_names = [
            i for i in self.ko_list_df.index + ": " + self.ko_list_df["gene_name"]
        ]
        self.show_genes = show_genes

        # Widgets
        self.filter_text = widgets.Text(
            value="",
            description="Filter KO list",
            placeholder="Type to filter..."
        )
        self.ko_formatted_names_w = widgets.Select(
            options=self.ko_formatted_names,
            description="KO list",
            rows=20,
            layout=widgets.Layout(width='600px')
        )
        self.ko_info_output = widgets.Output()

        # Widget behaviors
        self.filter_text.observe(self.filter_ko_list, names="value")
        self.ko_formatted_names_w.observe(self.display_ko_info, names="value")
        self.ko_formatted_names_w.value = None

        # Dashboard layout
        self.dashboard = widgets.HBox([
            widgets.VBox([self.filter_text, self.ko_formatted_names_w]),
            self.ko_info_output
        ])

    def validate_file_path(self, file_path, file_description):
        """Validate that the specified file path exists.

        Args:
            file_path (str): Path to the file.
            file_description (str): Description of the file being validated.

        Raises:
            FileNotFoundError: If the file does not exist.
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"{file_description} not found at path: {file_path}")

    def load_ko_list(self, ko_list_path, ko_synonyms_path):
        """Load and format the KO list with synonyms using Dask."""
        # Load KO list with Dask
        ko_list_df = dd.read_csv(ko_list_path, sep="\t")
        ko_list_df = ko_list_df.set_index(ko_list_df.columns[0])

        # Select the necessary column and rename it
        ko_list_df = ko_list_df[["definition"]]
        ko_list_df.columns = ["gene_name"]

        ko_synonyms = dd.read_csv(ko_synonyms_path, sep="\t", header=None, names=["synonyms"])
        ko_synonyms = ko_synonyms.fillna("")  # Dask equivalent of fillna

        # Merge the two Dask dataframes
        ko_list_df = ko_list_df.merge(ko_synonyms, left_index=True, right_index=True, how="left")
        ko_list_df["gene_name"] = ko_list_df["synonyms"] + ", " + ko_list_df["gene_name"]
        ko_list_df = ko_list_df.drop(columns=["synonyms"])

        # Compute the result and convert to pandas for further use
        ko_list_df = ko_list_df.compute()
        ko_list_df.index.name = "KO"
        return ko_list_df

    def filter_ko_list(self, change):
        """Filter KO list based on user input."""
        self.ko_formatted_names_w.value = None
        filter_value = self.filter_text.value.lower()
        filtered_options = [
            ko for ko in self.ko_formatted_names if filter_value in ko.lower()
        ]
        self.ko_formatted_names_w.options = filtered_options

        if not filter_value:
            self.ko_formatted_names_w.value = None

    def display_ko_info(self, change):
        """Fetch and display KO information for the selected KO."""
        self.ko_info_output.clear_output()
        selected_ko = change["new"]
        if selected_ko:
            ko_id = selected_ko.split(":")[0].strip()  # Extract KO ID
            with self.ko_info_output:
                # print(f"Fetching information for {ko_id}...")
                try:
                    ko_info = REST.kegg_get(f"ko:{ko_id}").read()
                    if self.show_genes==False:
                        ko_info = re.sub(r"GENES\s+.*?(?=\n\S)", "", ko_info, flags=re.DOTALL)
                    print(ko_info)
                except Exception as e:
                    print(f"Failed to fetch information for {ko_id}: {e}")

    def display(self):
        """Display the dashboard."""
        display(self.dashboard)