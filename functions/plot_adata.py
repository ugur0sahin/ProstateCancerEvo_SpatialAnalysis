import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import matplotlib as mpl


def plot_adata(adata_to_plot, cluster_key, plot_type="umap"):
    if plot_type == "spatial":

        coords = adata_to_plot.obsm['spatial'].copy()
        df = pd.DataFrame(coords, columns=['x', 'y'], index=adata_to_plot.obs_names)

        df[cluster_key] = adata_to_plot.obs[cluster_key].astype(str)

        tab20 = [mpl.colors.rgb2hex(c) for c in plt.get_cmap('tab20').colors]

        fig = px.scatter(
            df,
            x='x',
            y='y',
            color=cluster_key,
            title='Spatial scatter — manual cell types',
            #color_discrete_sequence=tab20,
            hover_name=df.index,
            width=1600,
            height=700
        )

        fig.update_traces(marker=dict(size=2, opacity=0.8))
        fig.update_yaxes(autorange='reversed')
        fig.update_layout(
            legend_title_text=cluster_key,
            legend=dict(
                itemsizing='constant',
                traceorder='normal',
                bgcolor='rgba(255,255,255,0.5)',
                x=1.02, y=1
            ),
            margin=dict(l=20, r=200, t=50, b=20)
        )

        fig.show()



    elif plot_type == "umap":

        df = pd.DataFrame(
            adata_to_plot.obsm['X_umap'],
            columns=['UMAP1', 'UMAP2'],
            index=adata_to_plot.obs_names
        ).copy()

        df[cluster_key] = adata_to_plot.obs[cluster_key].astype(str)

        tab20 = [mpl.colors.rgb2hex(c) for c in plt.get_cmap('tab20').colors]

        fig = px.scatter(
            df,
            x='UMAP1',
            y='UMAP2',
            color=cluster_key,
            title='UMAP embedding — Leiden clusters',
            color_discrete_sequence=tab20,
            hover_name=df.index,
            width=1400,
            height=1200
        )

        fig.update_traces(marker=dict(size=3, opacity=0.8))
        fig.update_layout(
            legend_title_text='Leiden cluster',
            legend=dict(
                itemsizing='constant',
                traceorder='normal',
                bgcolor='rgba(255,255,255,0.5)',
                x=1.02, y=1
            ),
            margin=dict(l=20, r=200, t=50, b=20)
        )

        fig.show()