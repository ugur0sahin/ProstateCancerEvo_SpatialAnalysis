import scanpy as sc
import spatialdata
import numpy as np
from scipy.stats import linregress
import pandas as pd


OUT = "/Users/ugursahin/ProstateCancerEvo_SpatialAnalysis/scripts/phenotype_associations/results"


for Tissue, TableFolder in [["Region1", "REGION1_TABLES"],["Region3", "REGION3_TABLES"],["Region4", "REGION4_TABLES"]]:
    
    adata_main = sc.read_h5ad(f"/Volumes/ProstateCancerEvoMain/dbs/Ongoing/{Tissue}/{TableFolder}/{Tissue}_Xenium_Phen_HE_Integrated.Xenium_Native_Table.Phen_and_ImagePropsEnriched.V2.h5ad")

    ## 
    adata_for_annotation = sc.read_h5ad(f"/Volumes/ProstateCancerEvoMain/dbs/Ongoing/{Tissue}/{TableFolder}/{Tissue}_Xenium_Phen_HE_Integrated.Xenium_Process_Table.Cell_Type_Annotated_Auto.V3.All_Genes.h5ad")
    adata_main.obs["cell_type_xenium_panel"] = adata_for_annotation.obs["cell_type_xenium_panel"].values


    for cell_type in set(adata_main.obs["cell_type_xenium_panel"].to_list()):

        mask = adata_main.obs["cell_type_xenium_panel"] == cell_type
        adata = adata_main[mask, :].copy()

        expr = adata.X

        if hasattr(expr, "toarray"):
            T = np.array(expr.toarray()).sum(axis=1)
        else:
            T = np.asarray(expr.sum(axis=1)).ravel()

        y = adata.obs["cell_area"].values


        records = []
        for goi in adata.var_names:
            
            x = np.asarray(adata[:, goi].X.toarray().flatten())
            if np.unique(x).size > 4:
                print(goi)
                # 2) define mask properly (discard zero-expression cells)
                mask = x > 0

                # 3) Option 3: use masked T and y
                T_y =  T / y
                x2, y2 = x[mask], T_y[mask]

                
                try:
                    slope, intercept, r_value, p_value, std_err = linregress(x2, y2)

                    # 95% CI for slope under normal approximation:
                    ci_low  = slope - 1.96 * std_err
                    ci_high = slope + 1.96 * std_err

                    # validation: CI does *not* cross zero
                    validation = (ci_low > 0) or (ci_high < 0)

                    records.append({
                    "goi":          goi,
                    "slope":        slope,
                    "intercept":    intercept,
                    "p_value":      p_value,
                    "r_value":      r_value,
                    "std_err":      std_err,
                    "ci_low":       ci_low,
                    "ci_high":      ci_high,
                    "validation":   validation
                    })
                except:
                    pass

        try:
            df = pd.DataFrame(records)
            print(df)
            df = df[df["validation"] == True]
            df.sort_values(by="slope")

            from statsmodels.stats.multitest import multipletests

            n_tests = len(df)
            reject, p_bonf, _, _ = multipletests(df['p_value'], alpha=0.05, method='bonferroni')

            df['p_bonf']        = p_bonf
            df['sig_bonferroni'] = reject

            print(f"{df['sig_bonferroni'].sum()} genes significant after Bonferroni correction")

            df[df['sig_bonferroni']].sort_values('slope')

            df = df[df["sig_bonferroni"] == True]
            print( f" For {Tissue} and {cell_type}, The # of genes has negative expression correlation with cell size {df[df['slope'] < 0].shape[0]}, and positive correlation  {df[df['slope'] > 0].shape[0]}")
            
            df.to_csv(f"{OUT}/{Tissue}__{cell_type}__FDR_CI_validated__Expression_CellArea_assoc.csv")
        except:
            print(f"Due to an erray {cell_type} in {Tissue} skipped ! ")