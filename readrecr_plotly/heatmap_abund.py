# just parts of the script, will require edits

import pandas as pd
import plotly.subplots as sp
import plotly.graph_objs as go
import numpy as np
import os
pd.set_option("display.max_columns", None)

reorder_idx = ["GCF_019343105_1", "metabat_73", "GCA_004298435_1", "GCF_004135935_1", "GCA_021780015_1", "GCA_021324915_1", "JOPF_NEAR_metabat2_504_refined", "JOPF_N14C_concoct_25_refined", "JOPF_N102_maxbin2_018_sub_refined", "JOPF_N143_metabat2_18_refined"]

reorder_cols = {
    "Bahram": ["ERR5766915", "ERR5766914", "ERR5766913", "ERR5766912", "ERR5766911", "ERR5766910", "ERR5766909", "ERR5766908", "ERR5766907", "ERR5766906", "ERR5766688", "ERR5766687", "ERR5766686", "ERR5766685", "ERR5766684", "ERR5766683", "ERR5766682", "ERR5766681", "ERR5766680", "ERR5766679", "ERR5766678", "ERR5766677", "ERR5766676", "ERR5766674", "ERR5766673", "ERR5766672", "ERR5766671", "ERR5766670", "ERR5766669", "ERR5766668", "ERR5766667", "ERR5766664", "ERR5766586", "ERR5766526", "ERR5766525", "ERR5766524", "ERR5766523", "ERR5766522", "ERR5766521", "ERR5766520", "ERR5766493", "ERR5766257", "ERR5763708", "ERR5763227", "ERR5763226", "ERR5763223", "ERR5763222", "ERR5762868", "ERR5762505", "ERR5762504", "ERR5762360", "ERR5762359", "ERR5762358", "ERR5762329", "ERR5762296", "ERR5762295", "ERR5762279", "ERR5762278", "ERR5762277", "ERR5762276", "ERR5762275", "ERR5762267", "ERR5762240", "ERR5762236", "ERR5762234", "ERR5762230", "ERR5762228", "ERR5762226", "ERR5762225", "ERR5762224", "ERR5762223", "ERR5762222", "ERR5762221", "ERR5762220", "ERR5762219", "ERR5762218", "ERR5762217", "ERR5762216", "ERR5762215", "ERR5762214", "ERR5762213", "ERR5762212", "ERR5762211", "ERR5762210", "ERR5762209", "ERR5762208", "ERR5762207", "ERR5762206", "ERR5762205", "ERR5762204", "ERR5762203", "ERR5762202", "ERR5762201", "ERR5762200", "ERR5762199", "ERR5761597", "ERR5761268", "ERR5758156", "ERR5758137", "ERR5758133", "ERR5758132", "ERR5757926", "ERR5757924", "ERR5757923", "ERR5757922", "ERR5757864", "ERR5757789", "ERR5757787", "ERR5757781", "ERR5757776", "ERR5757773", "ERR5757737", "ERR5757732", "ERR5757700", "ERR5757699", "ERR5757698", "ERR5757697", "ERR5757696", "ERR5757695", "ERR5757662", "ERR5757377", "ERR5757362", "ERR5757358", "ERR5757355", "ERR5757354", "ERR5756921", "ERR5756920", "ERR5756919", "ERR5754037", "ERR5754036", "ERR5754029", "ERR5754023", "ERR5754020", "ERR5754007", "ERR5754004", "ERR5753969", "ERR5753965", "ERR5753930", "ERR5753906", "ERR5753807", "ERR5753773", "ERR5753395", "ERR5753391", "ERR5753385", "ERR5753381", "ERR5753375", "ERR5752219", "ERR5752218", "ERR5752217", "ERR5751860", "ERR5751859", "ERR5751858", "ERR5751857", "ERR5751856", "ERR5750992", "ERR5750990", "ERR5750989", "ERR5750987", "ERR5750985", "ERR5750984", "ERR5750983", "ERR5750982", "ERR5750981", "ERR5750980", "ERR5750978", "ERR5750922", "ERR5750921", "ERR5750904", "ERR5750817", "ERR5750771", "ERR5750770", "ERR5750768", "ERR5750511", "ERR5750510", "ERR5749944", "ERR5749943", "ERR5749942", "ERR5749940", "ERR5749939", "ERR5749931", "ERR5749930", "ERR5749929", "ERR5749928", "ERR5749927", "ERR5749926", "ERR5749924", "ERR5749923", "ERR5749872", "ERR5749764", "ERR5748380", "ERR5748206", "ERR5748199", "ERR5748190", "ERR5748176", "ERR5748164"],
    "Bandla": ['SRR21679849', 'SRR21679850', 'SRR21679851', 'SRR21679852', 'SRR21679853', 'SRR21679854', 'SRR21679855', 'SRR21679856', 'SRR21679857', 'SRR21679858', 'SRR21679859', 'SRR21679860', 'SRR21679861', 'SRR21679862', 'SRR21679863', 'SRR21679864', 'SRR21679865', 'SRR21679866', 'SRR21679867', 'SRR21679868', 'SRR21679869', 'SRR21679870', 'SRR21679871', 'SRR21679872', 'SRR21679873', 'SRR21679874', 'SRR21679875', 'SRR21679876', 'SRR21679877', 'SRR21679878', 'SRR21679879', 'SRR21679880', 'SRR21679881', 'SRR21679882', 'SRR21679883', 'SRR21679884'],
    "Woodcroft": ['SRR23247659', 'SRR23247660', 'SRR23247661', 'SRR23247662', 'SRR23247663', 'SRR23247664', 'SRR23247665', 'SRR23247666', 'SRR23247667', 'SRR23247668', 'SRR23247669', 'SRR23247670', 'SRR23247671', 'SRR23247672', 'SRR23247673', 'SRR23247674', 'SRR23247675', 'SRR23247676', 'SRR23247677', 'SRR23247678', 'SRR23247679', 'SRR23247680', 'SRR23247681', 'SRR23247682', 'SRR23247683', 'SRR23247684', 'SRR23247685', 'SRR23247686', 'SRR23247687', 'SRR23247688', 'SRR23247689', 'SRR23247690', 'SRR23247691', 'SRR23247692', 'SRR23247693', 'SRR23247694', 'SRR23247695', 'SRR23247696', 'SRR23247697', 'SRR23247698', 'SRR23247699', 'SRR23247700', 'SRR23247701', 'SRR23247702', 'SRR23247703', 'SRR23247704', 'SRR23247705', 'SRR23247706', 'SRR23247707', 'SRR23247708', 'SRR23247709', 'SRR23247710', 'SRR23247711', 'SRR23247712', 'SRR23247713', 'SRR23247714', 'SRR23247715', 'SRR23247716', 'SRR23247717', 'SRR23247718', 'SRR23247719', 'SRR23247720', 'SRR23247721', 'SRR23247722', 'SRR23247723', 'SRR23247724', 'SRR23247725', 'SRR23247726', 'SRR23247727', 'SRR23247728', 'SRR23247729', 'SRR23247730', 'SRR23247731', 'SRR23247732', 'SRR23247733', 'SRR23247734', 'SRR23247735', 'SRR23247736', 'SRR23247737', 'SRR23247738', 'SRR23247739', 'SRR23247740', 'SRR23247741', 'SRR23247742', 'SRR23247743', 'SRR23247744', 'SRR23247745', 'SRR23247746', 'SRR23247747', 'SRR23247748', 'SRR23247749', 'SRR23247750', 'SRR23247751', 'SRR23247752', 'SRR23247753', 'SRR23247754', 'SRR23247755', 'SRR23247756', 'SRR23247757', 'SRR23247758', 'SRR23247759', 'SRR23247760', 'SRR23247761', 'SRR23247762', 'SRR23247763', 'SRR23247764', 'SRR23247765', 'SRR23247766', 'SRR23247767', 'SRR23247768', 'SRR23247769', 'SRR23247770', 'SRR23247771', 'SRR23247772', 'SRR23247773', 'SRR23247774', 'SRR23247775', 'SRR23247776', 'SRR23247777', 'SRR23247778', 'SRR23247779', 'SRR23247780', 'SRR23247781', 'SRR23247782', 'SRR23247783', 'SRR23247784', 'SRR23247785', 'SRR23247786', 'SRR23247787', 'SRR23247788', 'SRR23247789', 'SRR23247790', 'SRR23247791', 'SRR23247792', 'SRR23247793', 'SRR23247794', 'SRR23247795', 'SRR7151488', 'SRR7151489', 'SRR7151490', 'SRR7151491', 'SRR7151492', 'SRR7151493', 'SRR7151494', 'SRR7151495', 'SRR7151496', 'SRR7151497', 'SRR7151498', 'SRR7151499', 'SRR7151500', 'SRR7151501', 'SRR7151502', 'SRR7151503', 'SRR7151504', 'SRR7151505', 'SRR7151506', 'SRR7151507', 'SRR7151508', 'SRR7151509', 'SRR7151510', 'SRR7151511', 'SRR7151512', 'SRR7151513', 'SRR7151514', 'SRR7151515', 'SRR7151516', 'SRR7151517', 'SRR7151518', 'SRR7151519', 'SRR7151520', 'SRR7151521', 'SRR7151522', 'SRR7151523', 'SRR7151524', 'SRR7151525', 'SRR7151526', 'SRR7151527', 'SRR7151528', 'SRR7151529', 'SRR7151530', 'SRR7151531', 'SRR7151532', 'SRR7151533', 'SRR7151534', 'SRR7151535', 'SRR7151536', 'SRR7151537', 'SRR7151538', 'SRR7151539', 'SRR7151540', 'SRR7151541', 'SRR7151542', 'SRR7151543', 'SRR7151544', 'SRR7151545', 'SRR7151546', 'SRR7151547', 'SRR7151548', 'SRR7151549', 'SRR7151550', 'SRR7151551', 'SRR7151552', 'SRR7151553', 'SRR7151554', 'SRR7151555', 'SRR7151556', 'SRR7151557', 'SRR7151558', 'SRR7151559', 'SRR7151560', 'SRR7151561', 'SRR7151562', 'SRR7151563', 'SRR7151564', 'SRR7151565', 'SRR7151566', 'SRR7151567', 'SRR7151568', 'SRR7151569', 'SRR7151570', 'SRR7151571', 'SRR7151572', 'SRR7151573', 'SRR7151574', 'SRR7151575', 'SRR7151576', 'SRR7151577', 'SRR7151578', 'SRR7151579', 'SRR7151580', 'SRR7151581', 'SRR7151582', 'SRR7151583', 'SRR7151584', 'SRR7151585', 'SRR7151586', 'SRR7151587', 'SRR7151588', 'SRR7151589', 'SRR7151590', 'SRR7151591', 'SRR7151592', 'SRR7151593', 'SRR7151594', 'SRR7151595', 'SRR7151596', 'SRR7151597', 'SRR7151598', 'SRR7151599', 'SRR7151600', 'SRR7151601', 'SRR7151602', 'SRR7151603', 'SRR7151604', 'SRR7151605', 'SRR7151606', 'SRR7151607', 'SRR7151608', 'SRR7151609', 'SRR7151610', 'SRR7151611', 'SRR7151612', 'SRR7151613', 'SRR7151614', 'SRR7151615', 'SRR7151616', 'SRR7151617', 'SRR7151618', 'SRR7151619', 'SRR7151620', 'SRR7151621', 'SRR7151623', 'SRR7151624', 'SRR7151625', 'SRR7151626', 'SRR7151627', 'SRR7151628', 'SRR7151629', 'SRR7151630', 'SRR7151631', 'SRR7151632', 'SRR7151633', 'SRR7151634', 'SRR7151635', 'SRR7151636', 'SRR7151637', 'SRR7151638', 'SRR7151639', 'SRR7151640', 'SRR7151641', 'SRR7151642', 'SRR7151643', 'SRR7151644', 'SRR7151645', 'SRR7151646', 'SRR7151647', 'SRR7151648', 'SRR7151649', 'SRR7151650', 'SRR7151651', 'SRR7151652', 'SRR7151653', 'SRR7151654', 'SRR7151655', 'SRR7151656', 'SRR7151657', 'SRR7151658', 'SRR7151659', 'SRR7151660', 'SRR7151661', 'SRR7151662', 'SRR7151663', 'SRR7151664', 'SRR7151665', 'SRR7151666', 'SRR7151667', 'SRR7151668', 'SRR7151669', 'SRR7151670', 'SRR7151671', 'SRR7151672', 'SRR7151673', 'SRR7151674', 'SRR7151675', 'SRR7151676', 'SRR7151677', 'SRR7151678', 'SRR7151679', 'SRR7151680', 'SRR7151681', 'SRR7151682', 'SRR7151683', 'SRR7151684', 'SRR7151685', 'SRR7151686', 'SRR7151687', 'SRR7151688', 'SRR7151689', 'SRR7151690', 'SRR7151691', 'SRR7151692', 'SRR7151693', 'SRR7151694', 'SRR7151695', 'SRR7151696', 'SRR7151697', 'SRR7151698', 'SRR7151699', 'SRR7151700', 'SRR7151701', 'SRR7151702', 'SRR7151703', 'SRR7151704', 'SRR7151705', 'SRR7151706', 'SRR7151707', 'SRR7151708', 'SRR7151709', 'SRR7151710', 'SRR7151711', 'SRR7151712', 'SRR7151713', 'SRR7151714', 'SRR7151715', 'SRR7151716', 'SRR7151717', 'SRR7151718', 'SRR7151719', 'SRR7151720', 'SRR7151721', 'SRR7151722', 'SRR7151723', 'SRR7151724', 'SRR7151725', 'SRR7151726', 'SRR7151727', 'SRR7151728', 'SRR7151729', 'SRR7151730', 'SRR7151731', 'SRR7151732', 'SRR7151733', 'SRR7151734', 'SRR7151735', 'SRR7151736', 'SRR7151737', 'SRR7151738', 'SRR7151739', 'SRR7151740', 'SRR7151741', 'SRR7151742', 'SRR7151743'],
}

def get_dfm(mdf_path, abund_df, m):
    mdf = pd.read_csv(mdf_path, sep = "\t", index_col = 0)
    samples_in_df = [s for s in mdf.index.tolist() if s in abund_df.columns]
    
    dfm = pd.DataFrame(
        data = [mdf.loc[samples_in_df, m].tolist() for i in range(len(abund_df.index))],
            index = abund_df.index,
            columns = abund_df.columns
        )
    return dfm

# Load the abundance data
project = "Bahram"
abundance_df = pd.read_csv(f"anvio_summary/{project}/bins_across_samples/abundance.txt", sep="\t", index_col=0)
abundance_df = abundance_df.loc[reorder_idx, reorder_cols[f"{project}"]]

# Load the detection data
detection_df = pd.read_csv(f"anvio_summary/{project}/bins_across_samples/detection.txt", sep="\t", index_col=0)
detection_df = detection_df.loc[reorder_idx, reorder_cols[f"{project}"]]

mp = f"metadata/{project}_metadata.tsv"

hover_template = (
    "sample: %{x}<br>"
    "value: %{z}<br>"
    'site: %{customdata[0]}<br>'
    'land_use_1: %{customdata[1]}<br>'
    'landscape: %{customdata[2]}<br>'
    'vegetation: %{customdata[3]}<br>'
)

# Create customdata for both datasets
customdata_abundance = np.dstack((
    get_dfm(mp, abundance_df, "site").to_numpy(),
    get_dfm(mp, abundance_df, "land_use_1").to_numpy(),
    get_dfm(mp, abundance_df, "landscape").to_numpy(),
    get_dfm(mp, abundance_df, "vegetation").to_numpy(),
))

# Create a subplot with two rows
fig = sp.make_subplots(rows=1, cols=2, shared_yaxes=True)

# Add the abundance heatmap as the first trace
trace_abundance = go.Heatmap(
    z=abundance_df.values,
    x=abundance_df.columns,
    y=abundance_df.index,
    hovertemplate=hover_template,
    coloraxis = "coloraxis1",
    customdata=customdata_abundance,
)

# Add the detection heatmap as the second trace
trace_detection = go.Heatmap(
    z=detection_df.values,
    x=detection_df.columns,
    y=detection_df.index,
    hovertemplate=hover_template,
    coloraxis = "coloraxis2",
    customdata=customdata_abundance,  # Use the same customdata as abundance
    zmax = 1
)

# Add the traces to the subplots
fig.add_trace(trace_abundance, row=1, col=1)
fig.add_trace(trace_detection, row=1, col=2)

# Update the layout for both subplots
fig.update_layout(
    title=f"{project} Abundance and Detection Datasets",
    autosize=True,
    width=1600,
    height=800,  # Adjust the height as needed
    coloraxis=dict(colorscale='deep_r', colorbar_x=0.45, colorbar_thickness=23),
    coloraxis2=dict(colorscale='matter_r', colorbar_x=1, colorbar_thickness=23, cmin = 0, cmax = 1)
)

# Synchronize zooming behavior between the two subplots
fig.update_xaxes(matches='x')
fig.update_yaxes(matches='y')

# Customize the color bar length for each subplot
fig.update_coloraxes(colorbar_len=0.01, colorbar_yanchor='top', row=1, col=1)
fig.update_coloraxes(colorbar_len=0.01, colorbar_yanchor='top', row=1, col=2)

# Show the combined figure
fig.show()


###

def plot_everything(project):
    abundance_df = pd.read_csv(f"anvio_summary/{project}/bins_across_samples/abundance.txt", sep="\t", index_col=0)

    try:
        abundance_df = abundance_df.loc[reorder_idx, reorder_cols[f"{project}"]]
    except:
        print("Not reordering samples?")

    # Load the detection data
    detection_df = pd.read_csv(f"anvio_summary/{project}/bins_across_samples/detection.txt", sep="\t", index_col=0)

    try:
        detection_df = detection_df.loc[reorder_idx, reorder_cols[f"{project}"]]
    except:
        print("Not reordering samples?")

    # Create a subplot with two rows
    fig = sp.make_subplots(rows=1, cols=2, shared_yaxes=True)
    
    try:
        mp = f"metadata/{project}_metadata.tsv"

        hover_template = (
            "sample: %{x}<br>"
            "value: %{z}<br>"
            'cur_vegetation: %{customdata[0]}<br>'
            'Sample reference: %{customdata[1]}<br>'
        )

        # Create customdata for both datasets

        customdata_abundance = np.dstack((
            get_dfm(mp, abundance_df, "cur_vegetation").to_numpy(),
            get_dfm(mp, abundance_df, "Sample reference").to_numpy(),
        ))
        
        # Add the abundance heatmap as the first trace
        trace_abundance = go.Heatmap(
            z=abundance_df.values,
            x=abundance_df.columns,
            y=abundance_df.index,
            hovertemplate=hover_template,
            coloraxis="coloraxis1",
            customdata=customdata_abundance,
            zmin=0,  # Set the minimum value for the color scale
            zmax=1,  # Set the maximum value for the color scale
        )
        
        # Add the detection heatmap as the second trace
        trace_detection = go.Heatmap(
            z=detection_df.values,
            x=detection_df.columns,
            y=detection_df.index,
            hovertemplate=hover_template,
            coloraxis="coloraxis2",
            customdata=customdata_abundance,  # Use the same customdata as abundance
            zmin=0,  # Set the minimum value for the color scale
            zmax=1,  # Set the maximum value for the color scale to 1
        )
    except:
        print("Metadata missing")
            # Add the abundance heatmap as the first trace
        a_hover_template = (
            "sample: %{x}<br>"
            "abundance: %{z}<br>"
        )
        # Add the abundance heatmap as the first trace
        trace_abundance = go.Heatmap(
            z=abundance_df.values,
            x=abundance_df.columns,
            y=abundance_df.index,
            hovertemplate=a_hover_template,
            coloraxis="coloraxis1",
            zmin=0,  # Set the minimum value for the color scale
            zmax=1,  # Set the maximum value for the color scale
        )
        
        d_hover_template = (
            "sample: %{x}<br>"
            "detection: %{z}<br>"
        )
        # Add the detection heatmap as the second trace
        trace_detection = go.Heatmap(
            z=detection_df.values,
            x=detection_df.columns,
            y=detection_df.index,
            hovertemplate=d_hover_template,
            coloraxis="coloraxis2",
            zmin=0,  # Set the minimum value for the color scale
            zmax=1,  # Set the maximum value for the color scale to 1
        )


    # Add the traces to the subplots
    fig.add_trace(trace_abundance, row=1, col=1)
    fig.add_trace(trace_detection, row=1, col=2)

    # Update the layout for both subplots
    fig.update_layout(
        title=f"{project} Abundance and Detection Datasets",
        autosize=True,
        width=1600,
        height=800,  # Adjust the height as needed
        coloraxis=dict(colorscale='deep_r', colorbar_x=0.45, colorbar_thickness=23),
        coloraxis2=dict(colorscale='matter_r', colorbar_x=1, colorbar_thickness=23, cmin = 0, cmax = 1)
    )

    # Synchronize zooming behavior between the two subplots
    fig.update_xaxes(matches='x')
    fig.update_yaxes(matches='y')

    # Customize the color bar length for each subplot
    fig.update_coloraxes(colorbar_len=0.01, colorbar_yanchor='top', row=1, col=1)
    fig.update_coloraxes(colorbar_len=0.01, colorbar_yanchor='top', row=1, col=2)

    # Show the combined figure
    fig.show()