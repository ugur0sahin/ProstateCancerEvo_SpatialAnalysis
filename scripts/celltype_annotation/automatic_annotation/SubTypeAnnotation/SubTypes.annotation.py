import pandas as pd
import scanpy as sc
import spatialdata
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import decoupler as dc

Majors = ["ImmuneAnn_SubTypes","Mesenchymal_Subtypes","Epithelia_SubTypes"]

for Tissue in ["Region1","Region2","Region3","Region4"]:
    MajorSets_P = []
    for MajorSet in Majors:
        xen_file = f"/Volumes/ProstateCancerEvoMain/dbs/Ongoing/Collection/MajorAnnGroups/{MajorSet}/{Tissue}.{MajorSet}.h5ad"
        MajorSets_P.append(xen_file.replace(".h5ad",".AnnDC.h5ad"))
        adata_xen = sc.read_h5ad(xen_file) 
        adata_xen

        ann_file = f"/Volumes/ProstateCancerEvoMain/dbs/Ongoing/Collection/MajorAnnGroups/{MajorSet}/{MajorSet}.V1.xlsx"
        outpath = "/Volumes/ProstateCancerEvoMain/dbs/Ongoing/Collection/MajorAnnGroups/Results"
        df = pd.read_excel(ann_file)



        df["gene_list"] = df["geneSymbolmore1"].str.split(",")
        df_long = df.explode("gene_list")
        df_long["gene_list"] = df_long["gene_list"].str.strip()
        df_long = df_long[["cellName", "gene_list"]].rename(columns={"gene_list": "geneSymbol"})
        markers = df_long.reset_index(drop=True)

        

        dup_mask = markers["geneSymbol"].duplicated(keep=False)

        markers = markers[~dup_mask].copy()



        dc.run_ora(
            mat=adata_xen,
            net=markers,
            source='cellName',
            target='geneSymbol',
            min_n=3,
            verbose=True,
            use_raw=False
        )

        adata_xen.obsm['ora_estimate']
        acts = dc.get_acts(adata_xen, obsm_key='ora_estimate')

        # We need to remove inf and set them to the maximum value observed for pvals=0
        acts_v = acts.X.ravel()
        max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
        acts.X[~np.isfinite(acts.X)] = max_e

        res = dc.rank_sources_groups(acts, groupby='leiden', reference='rest', method='t-test_overestim_var')

        n_ctypes = 3
        ctypes_dict = res.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
        ctypes_dict

        # If adata_xen.obs['leiden'] is integer-typed, convert to string so it matches ctypes_dict keys:
        if pd.api.types.is_integer_dtype(adata_xen.obs['leiden']):
            adata_xen.obs['leiden'] = adata_xen.obs['leiden'].astype(str)

        # 4) Create a new column "assigned_celltype" by mapping each cell's cluster to ctypes_dict[cluster][0]
        def pick_first_celltype(cluster_label):
            # Return the first element of the list in ctypes_dict, if the key exists
            if cluster_label in ctypes_dict:
                return ctypes_dict[cluster_label][0]
            else:
                return None

        adata_xen.obs['assigned_celltype'] = adata_xen.obs['leiden'].map(pick_first_celltype)



        adata_xen.write_h5ad(xen_file.replace(".h5ad",".AnnDC.h5ad"))
        adata_xen.write_h5ad(
        f"{outpath}/{Tissue}/{xen_file.split('/')[-1].replace('.h5ad','.AnnDC.h5ad')}")


    adata1, adata2, adata3 = sc.read_h5ad(MajorSets_P[0]), sc.read_h5ad(MajorSets_P[1]), sc.read_h5ad(MajorSets_P[2])
    
    adata_concat = sc.concat(
    [adata1, adata2, adata3],
    join="outer",             # "inner" for intersection, "outer" for union of var_names 
    label="MajorTypes",            # name of the new obs‐column that stores which adata each cell came from
    keys=Majors,   
                              # these keys will populate: adata_concat.obs["batch"] = "batch1"/"batch2"/"batch3"
    index_unique="_",         # if you need to make obs_names unique, a suffix "_" will be appended
    )
    
    adata_concat.write_h5ad(f"/Volumes/ProstateCancerEvoMain/dbs/Ongoing/Collection/MajorAnnGroups/Results/{Tissue}/{Tissue}.AnnDc.All.h5ad")