import os
import pandas as pd
import ipywidgets as widgets
import plotly.express as px
from IPython.display import clear_output

ko_list_df_path = "/home/zarul/Zarul/Peatlands/status_nanopore/Methylovirgula/dram_faa/ko_list.txt"
ko_list_df = pd.read_csv(ko_list_df_path, sep="\t", header=None, index_col=0)
ko_list_df.columns = ["gene_name"]
ko_list_df.index.name = "KO"
ko_list = (ko_list_df.index + ": " + ko_list_df["gene_name"]).to_dict()

ko_formatted_names = (ko_list_df.index + ": " + ko_list_df["gene_name"]).tolist()
max_char = 50
ko_formatted_names_shortened = []

for i in ko_formatted_names:
    if len(i) > max_char:
        ko_formatted_names_shortened.append(f"{i[:max_char]}..")
    else:
        ko_formatted_names_shortened.append(i)
    
ko_formatted_names_w = widgets.Select(
    options=ko_formatted_names,
    description="KO list",
    rows=40,
    layout=widgets.Layout(width='450px')
)

def read_df(sample):
    df = pd.read_csv(f"dram_faa/{sample}/annotations.tsv", sep="\t")
    return df

def count_sample_ko(annot_df, ko):
    return annot_df.loc[annot_df["ko_id"] == ko, :].shape[0]

def matrix_ko(genome_selection, gene_selection):
    df_count_ko = pd.DataFrame(index=genome_selection)

    for gene in gene_selection:
        gene_name = ko_list[gene]
        if len(gene_name) > max_char:
            gene_name = f"{gene_name[:max_char]}..."
        counts = {}
        for genome in genomes_w.value:
            df = read_df(genome)
            counts[genome] = count_sample_ko(df, gene)

        df_count_ko[gene_name] = df_count_ko.index.map(counts)

    return df_count_ko

with open("genomes.txt", "r") as f:
    genomes = f.read()
genomes = genomes.split()

genomes_w = widgets.SelectMultiple(
    options=genomes,
    value=genomes,
    rows=40,
    description='Genomes',
    disabled=False,
    layout=widgets.Layout(max_width='400px')
)

init_ko_list = ["K00370", "K00371", "K00374", "K00373"]
ko_input = widgets.TagsInput(
    value=init_ko_list,
    description="Give some KEGG KOs",
    allow_duplicates=False
)

filter_text = widgets.Text(
    value="",
    description="Filter KO list",
    placeholder="Type to filter..."
)

# Create an interactive output for the heatmap
heatmap_output = widgets.Output()

colors = [
    "#e8e8e8",
    "#df0024",
    "#0085c7",
]

# Function to update the heatmap based on widget values
def update_heatmap(change):
    with heatmap_output:
        # Clear previous output
        clear_output(wait=True)

        # Get the selected genomes and KOs
        selected_genomes = genomes_w.value
        selected_kos = ko_input.value

        # Generate the matrix based on selected genomes and KOs
        mk = matrix_ko(selected_genomes, selected_kos)
        mk_max = mk.max().max()
        mk_val1 = 1/mk_max
        # Plot the heatmap
        fig = px.imshow(mk, aspect="auto")
        fig.update_layout(height=800, width=1400)
        fig.update_xaxes(tickangle=-45)

        if mk_val1 < 1:
            colorscale=[
                    (0.0, colors[0]),
                    (mk_val1-(0.25*mk_val1), colors[0]),
                    (mk_val1-(0.25*mk_val1), colors[1]),
                    (mk_val1+(0.25*mk_val1), colors[1]),
                    (mk_val1+(0.25*mk_val1), colors[2]),
                    (1, colors[2]),
                ]
        else:
            colorscale=[
                (0.0, colors[0]),
                (0.75, colors[0]),
                (0.75, colors[1]),
                (1, colors[1]),
            ]

        fig.update_coloraxes(
            colorscale=colorscale,
        )
        fig.show(config={"editable": True})


# Register the update function with the widgets
genomes_w.observe(update_heatmap, names='value')
ko_input.observe(update_heatmap, names='value')

# Function to filter KO list based on text input
def filter_ko_list(change):
    # Get the filter text
    filter_value = filter_text.value.lower()
    ko_formatted_names_w.value = None
    # Filter KO list based on the text input
    filtered_ko_list = []
    for ko in ko_formatted_names:
        if filter_value in ko.lower():
            filtered_ko_list.append(ko)

    # Update the KO formatted names widget
    ko_formatted_names_w.options = filtered_ko_list

# Register the filter function with the text widget
filter_text.observe(filter_ko_list, names='value')

# Function to handle copying selected item to ko_input widget
def add_to_ko_input(change):
    selected_item = ko_formatted_names_w.value
    if selected_item != None:
        if ko_input.value != None:
            if selected_item not in ko_input.value:
                mod_selected_item = selected_item.split(":")[0]
                ko_input.value = ko_input.value + [mod_selected_item]
        # print(f"Item '{selected_item}' added to ko_input.")

# Register the add function with the widget
ko_formatted_names_w.observe(add_to_ko_input, names='value')

# Display the initial heatmap
update_heatmap(None)


ko_input_label = widgets.HTML(
    value="<b>Selected KEGG KO</b>",
)


# Display the widgets and the interactive output
widgets.VBox(
    [widgets.HBox(
        [genomes_w,    heatmap_output,    widgets.VBox([filter_text, ko_formatted_names_w])]
    ),
    ko_input_label,
    ko_input]
)