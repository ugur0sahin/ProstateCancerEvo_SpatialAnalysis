import sopa
import spatialdata
import sopa.aggregation as agg
import pandas as pd
import scanpy as sc
import numpy as np
from skimage import data, filters, measure, morphology
import plotly.express as px
from skimage.draw import polygon
import plotly.graph_objects as go
import tifffile as tfi
import os
import shutil
import scipy.sparse as sp
from collections import Counter

""" ####################################################### """
""" _____________  Introducing The Parameters ______________"""
""" ####################################################### """


Tissue = "Region1"
TableFolder = "REGION1_TABLES"

out_path = "/Volumes/ProstateCancerEvo_SpatialAnalysis/dbs/Ongoing/{Tissue}/{TableFolder}"
os.makedirs(out_path, exist_ok=True)

channel_names = {
    '0':  'DAPI_PhenCyc', 
    '1':  'CD8_PhenCyc', 
    '2':  'CD31_PhenCyc', 
    '3':  'CD11c_PhenCyc', 
    '4':  'CD34_PhenCyc', 
    '5':  'CD20_PhenCyc', 
    '6':  'CD4_PhenCyc', 
    '7':  'CD21_PhenCyc', 
    '8':  'Ki67_PhenCyc', 
    '9':  'CD68_PhenCyc', 
    '10': 'E_cadherin_PhenCyc', 
    '11': 'pten_PhenCyc', 
    '12': 'NF-H_PhenCyc', 
    '13': 'FOXp3_PhenCyc', 
    '14': 'CD45_PhenCyc', 
    '15': 'TFAM_PhenCyc', 
    '16': 'PDL-1_PhenCyc', 
    '17': 'SOX2_PhenCyc', 
    '18': 'Vimentin_PhenCyc', 
    '19': 'HIF1A_PhenCyc', 
    '20': 'CD44_PhenCyc', 
    '21': 'Keratin_8_18_PhenCyc', 
    '22': 'aSMA_PhenCyc', 
    '23': 'EpCAM_PhenCyc', 
    '24': 'K14_PhenCyc'
}

proportion_constant = 4.7065
Xenium_Image = "Xenium_Image"

image_names = {

    "HE" : "HE_Image",
    #"morphology_focus" : "Xenium_Image",
    #"postXenium_5k_H03.ome" : "PhenoCycler_Image"
}

default_table_name = "table"
default_cell_boundaries_name = "cell_boundaries"

deleted_elements = []
#deleted_elements = ["Xenium_Lnmks", "XeniumMrks"]


required_columns = ['label',
                     'area', 
                     "centroid-0",
                     "centroid-1",
                     "centroid_local-1",
                     "centroid_weighted_local-0", "centroid_weighted_local-1"  ,
                     "euler_number",
                     "extent",
                     "feret_diameter_max",
                     "intensity_max",
                     "intensity_min",
                     "intensity_mean",
                     "solidity",
                     'eccentricity',
                     'perimeter',
                     "equivalent_diameter_area"]


path = "/Volumes/ProstateCancerEvo_SpatialAnalysis/dbs/Region_1_Xenium_Phen_HE_Integrated.V1.zarr"
sdata = spatialdata.read_zarr(path)


""" ####################################################### """
""" __________ Some Element Name Manipulation ______________"""
""" ####################################################### """


for k, v in image_names.items():
    sdata.images[v] = sdata.images[k]
print(sdata)


sdata.tables["GeneTranscripts_XStock_Native"] = sdata.tables[default_table_name]
sdata.tables["GeneTranscripts_XStockP"] = sdata.tables[default_table_name]
sdata.tables["GeneTranscripts_XStock_Backup"] = sdata.tables[default_table_name]

sdata.shapes["cell_boundaries_Backup"] = sdata.shapes[default_cell_boundaries_name]


for k, v in image_names.items():
    del sdata.images[k]
    
if deleted_elements != []:
    for el in deleted_elements:
        del sdata.points[el]

del sdata.tables[default_table_name]

""" ####################################################### """
""" ______Starting Point of Aggregation Process ____________"""
""" ####################################################### """

agg.aggregate(
    sdata,
    image_key="Xenium_Image",
    aggregate_genes=True,              
    points_key="transcripts",   
    shapes_key="cell_boundaries_Backup",
    min_intensity_ratio = 0,
    key_added="GeneTranscripts_XStockC"
)

agg.aggregate(
    sdata,
    aggregate_channels=True,
    aggregate_genes = False,
    image_key="PhenoCycler_Image",
    shapes_key="cell_boundaries_Backup",
    min_intensity_ratio = 0,
    key_added="Proteins_PhenoCycler_byXStock"              
)


""" ########################################################### """
""" ______ Renaming The Index and Cell_IDs w/ Xenium Stock _____"""
""" ########################################################### """

adata_stock = sdata.tables["GeneTranscripts_XStock_Native"].copy()

for table in ["Proteins_PhenoCycler_byXStock", "GeneTranscripts_XStockC"]:
    adata_ob = sdata.tables[table].copy()
    adata_ob.obs["cell_id"] = adata_stock.obs["cell_id"].to_list()
    adata_ob.index = adata_stock.obs["cell_id"].index

    for feat in ["region", "slide", "area"]:
        adata_ob.obs[feat].index = adata_stock.obs["cell_id"].to_list()
        
    sdata.tables[table] = adata_ob


""" ############################################################## """
""" _____________________ Record Them All _________________________"""
""" ############################################################## """


sdata.write(f"/Volumes/ProstateCancerEvo_SpatialAnalysis/dbs/Ongoing/{Tissue}/{Tissue}_Xenium_Phen_HE_Integrated.V2.ongoing.zarr")

adata_phen = sdata.tables["Proteins_PhenoCycler_byXStock"]
adata_phen.write_h5ad(f"{out_path}/{Tissue}_Xenium_Phen_HE_Integrated.Protein_PhenCycTable.V1.h5ad")

adata_native = sdata.tables["GeneTranscripts_XStock_Native"]
adata_native.write_h5ad(f"{out_path}/{Tissue}_Xenium_Phen_HE_Integrated.GeneTranscripts_XStock_Native.V1.h5ad")

adata_xenium_p = sdata.tables["GeneTranscripts_XStockP"]
adata_xenium_p.write_h5ad(f"{out_path}/{Tissue}_Xenium_Phen_HE_Integrated.GeneTranscripts_XStock_Process.V1.h5ad")

sdata.tables["GeneTranscripts_XStock_Backup"].write_h5ad(f"{out_path}/{Tissue}_Xenium_Phen_HE_Integrated.GeneTranscripts_XStock_Backup.V1.h5ad")
sdata.tables["GeneTranscripts_XStockC"].write_h5ad(f"{out_path}/{Tissue}_Xenium_Phen_HE_Integrated.GeneTranscripts_XStock_Calculated.V1.h5ad")


""" ######################################################################################## """
""" _____________________ Image Feature Extraction from Cell Shapes _________________________"""
""" ######################################################################################## """

## ReGenerate Cell Shapes in terms of Polygons !
new_ids = adata_native.obs["cell_id"].tolist()
shape = sdata.shapes['cell_boundaries']

records = []
for new_id, (old_id, poly) in zip(new_ids, shape['geometry'].items()):
    for vi, (x, y) in enumerate(poly.exterior.coords, start=1):
        records.append({
            'cell_id': new_id,
            'vertex':  vi,
            'x':       x,
            'y':       y
        })

shape_df = pd.DataFrame(records)


## ReGenerate Cell Shapes in terms of Polygons !
C = 0 # Don't need to change just used to calculate width x height of matrix
img = sdata.images[Xenium_Image]['scale0']['image'].isel(c=C)
arr    = np.asarray(img)
H, W = arr.shape
mask = np.zeros((H, W), dtype=np.int32)


labels = {}
for label, (cell_id, grp) in enumerate(shape_df.groupby('cell_id'), start=1):
    labels[label] = cell_id
    cols = np.floor(grp['x'].values * proportion_constant).astype(int) # Proportion constant must be adjusted due to conversion from pixel to global
    rows = np.floor(grp['y'].values * proportion_constant).astype(int)

    cols = np.clip(cols, 0, W-1)
    rows = np.clip(rows, 0, H-1)

    rr, cc = polygon(rows, cols, shape=mask.shape)
    mask[rr, cc] = label

os.makedirs(out_path, exist_ok=True)
np.savetxt(f"tmp/{Tissue}_mask.csv", mask, fmt="%d", delimiter=",")
records_df = pd.DataFrame(records)
records_df.to_csv(f"tmp/{Tissue}_mask_labels.csv", header=False, index=False)


# Start for Features We Interest
for i in range(0, 4):
    try:
        img_raw = sdata.images[Xenium_Image]['scale0']['image'].isel(c=i)
        label_image = mask
        
        img = np.asarray(img_raw)

        props = measure.regionprops_table(
            label_image,
            intensity_image=img,
            properties=['area',
        'area_bbox',
        'area_convex',
        'area_filled',
        'axis_major_length',
        'axis_minor_length',
        'centroid',
        'centroid_local',
        'centroid_weighted',
        'centroid_weighted_local',
        'coords',
        'eccentricity',
        'equivalent_diameter_area',
        'euler_number',
        'extent',
        'feret_diameter_max',
        'image',
        'image_convex',
        'image_filled',
        'image_intensity',
        'inertia_tensor',
        'inertia_tensor_eigvals',
        'intensity_max',
        'intensity_mean',
        'intensity_min',
        'label',
        'moments',
        'moments_central',
        'moments_hu',
        'moments_normalized',
        'moments_weighted',
        'moments_weighted_central',
        'moments_weighted_hu',
        'moments_weighted_normalized',
        'orientation',
        'perimeter',
        'perimeter_crofton',
        'slice',
        'solidity']
        )


        props_df = pd.DataFrame(props)
        props_df.to_csv(f"{Tissue}_ExtractedProps"+"_forChannel"+str(i)+".csv", header=True)
    except:
        print("Error in processing image for channel: ", i)
        continue



""" ######################################################################################## """
""" ___________ Attribution of Features to ADATA Native to undergo Shape ~ Gene  ____________"""
""" ######################################################################################## """



adata_b = adata_native.copy()
adata_a = adata_phen.copy()


# 2) Extract a dense array from adata_a.X if necessary
X = adata_a.X.toarray() if sp.issparse(adata_a.X) else adata_a.X

# 4) Loop over each column and assign it using the channel_names dict
for idx in range(X.shape[1]):
    col_name = channel_names[str(idx)]
    adata_b.obs[col_name] = X[:, idx]

adata = adata_b


#Read all of the properties CSVs and collect into one merged DataFrame

# 1) Read all your CSVs into a list of DataFrames
dfs = []
for i in range(0,2):
    path = f"{Tissue}_ExtractedProps_forChannel{i}.csv"
    # assume your files have a header row and the first column is an index
    df = pd.read_csv(path, index_col=0)
    dfs.append(df)

# 2) Make sure they all have the same shape, index & columns
assert all(df.shape == dfs[0].shape for df in dfs)
assert all((df.index == dfs[0].index).all()   for df in dfs)
assert all((df.columns == dfs[0].columns).all() for df in dfs)

# 3) Create the merged DataFrame of object‐dtype
merged = dfs[0].copy().astype(object)

# 4) Fill each cell with a list of values from each df
for idx in merged.index:
    for col in merged.columns:
        merged.at[idx, col] = [df.at[idx, col] for df in dfs]


# 4b) After you’ve built each list, check for a majority
for idx in merged.index:
    for col in merged.columns:
        vals = [df.at[idx, col] for df in dfs]  # your list of 4 values
        cnt = Counter(vals)
        # find the most common value and its count
        most_common_val, count = cnt.most_common(1)[0]
        if count > len(vals) / 2:
            # majority found → replace with the single value
            merged.at[idx, col] = most_common_val
        else:
            # no majority → keep the full list
            merged.at[idx, col] = vals


props_df_req = merged[required_columns].iloc[:, 1:]
props_df_req.index = props_df_req.index.astype(str)
obs_props_df_req = adata.obs.join(props_df_req)
obs_props_df_req


adata.obs = obs_props_df_req

# 1. Find all object-dtype columns in obs
obj_cols = [c for c, dt in adata.obs.dtypes.items() if dt == "object"]
print("Object-dtype columns:", obj_cols)

# 2. For each such column, inspect what types you actually have
for c in obj_cols:
    sample_types = set(type(v) for v in adata.obs[c].iloc[:100])
    print(f"  {c}: {sample_types}")


# 3b. Or, if you want to keep them, convert to string:
adata.obs[obj_cols] = adata.obs[obj_cols].astype(str)

# 4. (Optional) ensure index and column names are strings
adata.obs.index   = adata.obs.index.map(str)
adata.obs.columns = adata.obs.columns.map(str)

adata.write_h5ad(f"{out_path}/{Tissue}_Xenium_Phen_HE_Integrated.Xenium_Native_Table.Phen_and_ImagePropsEnriched.V2.h5ad")
shutil.rmtree("tmp/{Tissue}_mask_labels.csv", ignore_errors=True)
shutil.rmtree("tmp/{Tissue}_mask.csv", ignore_errors=True)