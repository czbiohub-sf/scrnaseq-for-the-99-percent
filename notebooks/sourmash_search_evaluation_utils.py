
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import roc_curve, roc_auc_score




def pivot_predicted_similarity_scores(
    search_results, 
    celltype_col, 
    sbt_organism,
    query_organism,
    scoring_groupby,
):
    predicted_celltype_col = f'{sbt_organism}_{celltype_col}'
    groundtruth_celltype_col = f'{query_organism}_{celltype_col}'

    similarity_scores = search_results.pivot_table(
        index=[f'{query_organism}_cell_id'],
        columns=[predicted_celltype_col] + scoring_groupby,
        values='similarity', aggfunc=np.max
    )
    return similarity_scores

def get_celltype_dummies(
    adata, 
    celltype_col, 
    sbt_organism,
    query_organism,
):
    predicted_celltype_col = f'{sbt_organism}_{celltype_col}'
    groundtruth_celltype_col = f'{query_organism}_{celltype_col}'
    
    # remove 1 from e.g. "mouse1" --> "mouse"
    query_organism_no_number = query_organism.strip('1234567890')

    ground_truth = pd.get_dummies(adata.obs.query('species_batch == @query_organism_no_number')[celltype_col])
    return ground_truth


def compute_roc_per_celltype(
    pivoted_similarity_scores, celltype_dummies,
        celltype_col, 
    sbt_organism,
    query_organism,
):
    
    dfs = []

    for compartment, y_pred in pivoted_similarity_scores.iteritems():
        y_true = celltype_dummies.loc[y_pred.index, compartment]
        fpr, tpr, thresholds = roc_auc_score(y_true, y_pred)
        df = pd.DataFrame({'fpr': fpr, 'tpr': tpr, 'thresholds': thresholds})
        df[celltype_col] = compartment
        dfs.append(df)
    roc_df = pd.concat(dfs, ignore_index=True)
    return roc_df


def compute_roc_auc_per_celltype(
    pivoted_similarity_scores, celltype_dummies,
        celltype_col, 
    sbt_organism,
    query_organism,
):
    
    dfs = []

    for compartment, y_pred in pivoted_similarity_scores.iteritems():
        y_true = celltype_dummies.loc[y_pred.index, compartment]
        fpr, tpr, thresholds = roc_(y_true, y_pred)
        df = pd.DataFrame({'fpr': fpr, 'tpr': tpr, 'thresholds': thresholds})
        df[celltype_col] = compartment
        dfs.append(df)
    roc_df = pd.concat(dfs, ignore_index=True)
    return roc_df

def plot_roc_per_celltype(roc_df, celltype_col, **kwargs):

    g = sns.FacetGrid(data=roc_df, hue=celltype_col, aspect=1.1, **kwargs)
    g.map(plt.plot, 'fpr', 'tpr')
    g.add_legend()
    return g


def compute_and_plot_roc(
    search_results, celltype_col, 
    sbt_organism,
    query_organism, 
    adata,
    palette='tab20', **kwargs):

    pivoted_similarity = pivot_predicted_similarity_scores(
        search_results, 
                celltype_col=celltype_col, 
        sbt_organism=sbt_organism,
        query_organism=query_organism,
    )
    celltype_dummies = get_celltype_dummies(
        adata,
                celltype_col=celltype_col, 
        sbt_organism=sbt_organism,
        query_organism=query_organism,
    )
    roc_df = compute_roc_per_celltype(
        pivoted_similarity, celltype_dummies,
                celltype_col=celltype_col, 
        sbt_organism=sbt_organism,
        query_organism=query_organism,
    )
    g = plot_roc_per_celltype(
        roc_df,         
                          celltype_col=celltype_col, palette=palette, **kwargs)
    return g
