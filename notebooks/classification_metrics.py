from sklearn.metrics import (
    adjusted_rand_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
)
import pandas as pd


def dataframeize(series, score_name):
    df = series.reset_index()
    df = df.rename(columns={0: "score_value"})
    df["score_name"] = score_name
    return df


def get_adjusted_rand_scores(
    predicted_cells,
    scoring_groupby,
    ground_truth_celltype_col,
    predicted_celltype_col,
    # Take labels as an input, but don't do anything with it
    # This is to be comparible with the inputs to get_f1_scores without adding special cases
    labels=None
):
    adjusted_rand_scores = predicted_cells.groupby(scoring_groupby).apply(
        lambda x: adjusted_rand_score(
            x[ground_truth_celltype_col],
            x[predicted_celltype_col],
        )
    )
    adjusted_rand_scores_df = dataframeize(adjusted_rand_scores, "adjusted_rand_score")
    return adjusted_rand_scores_df


def get_f1_scores(
    predicted_cells,
    scoring_groupby,
    ground_truth_celltype_col,
    predicted_celltype_col,
    labels=None
):
    f1_scores = predicted_cells.groupby(scoring_groupby).apply(
        lambda x: f1_score(
            x[ground_truth_celltype_col],
            x[predicted_celltype_col],
            sample_weight=x.similarity,
            average="weighted",
            labels=labels
        )
    )
    f1_scores_df = dataframeize(f1_scores, "f1_score")
    return f1_scores_df


scorers = (get_f1_scores, get_adjusted_rand_scores)


def get_f1_ari_scores(
    predicted_cells,
    ground_truth_celltype_col,
    predicted_celltype_col,
    scoring_groupby,
    scorers=scorers,
):

    dfs = []
    for scorer in scorers:
        df = scorer(
            predicted_cells,
            scoring_groupby=scoring_groupby,
            ground_truth_celltype_col=ground_truth_celltype_col,
            predicted_celltype_col=predicted_celltype_col,
        )
        dfs.append(df)
    concated_scores_df = pd.concat(dfs, ignore_index=True)
    return concated_scores_df


def get_ksize_maximizing_mean(classification_metrics):
    score_value_means = classification_metrics.groupby(
        ["alphabet", "score_name", "ksize"]
    ).score_value.mean()
    score_value_means_ksize_argmax = score_value_means.groupby(
        level=[0, 1], group_keys=False
    ).apply(lambda x: x.nlargest(1))
    score_value_means_ksize_argmax = score_value_means_ksize_argmax.reset_index()
    score_value_means_ksize_argmax[
        "mean_score"
    ] = score_value_means_ksize_argmax.score_name.str.split("_score").str[0]
    score_value_means_ksize_argmax[
        "alphabet_ksize"
    ] = score_value_means_ksize_argmax.apply(
        lambda x: f"{x.alphabet}, ksize: {x.ksize}", axis=1
    )
    return score_value_means_ksize_argmax


def pivot_and_plot_argmax_ksize_scores(score_value_means_ksize_argmax):
    score_value_means_ksize_argmax_2d = score_value_means_ksize_argmax.pivot(
        index="alphabet_ksize", columns="mean_score", values="score_value"
    )
    fig, ax = plt.subplots(figsize=(3, 2))
    sns.heatmap(
        score_value_means_ksize_argmax_2d,
        annot=True,
        vmin=0.5,
        vmax=1,
        cmap="viridis_r",
        fmt=".5g",
    )
    ax.set(ylabel=None)


def subsample_and_score(
    predicted_cells_subset,
    predicted_celltype_col,
    ground_truth_celltype_col,
    scoring_groupby,
    n_iterations=1000,
    n_cells_per_sample=3,
    labels=None,
    scorers=(get_f1_scores, get_adjusted_rand_scores),
    keys=None,
):
    score_dfs = []
    for i in range(n_iterations):
        subsampled = predicted_cells_subset.groupby(predicted_celltype_col).apply(
            lambda x: x.sample(n_cells_per_sample)
            if len(x) >= n_cells_per_sample
            else None
        )
        for scorer in scorers:
            score_df = scorer(
                subsampled,
                predicted_celltype_col=predicted_celltype_col,
                ground_truth_celltype_col=ground_truth_celltype_col,
                scoring_groupby=scoring_groupby,
                labels=labels
            )
            score_df["iteration"] = i

            for name, key in zip(scoring_groupby, keys):
                score_df[name] = key

            score_dfs.append(score_df)

    concat_score_df = pd.concat(score_dfs, ignore_index=True)
    return concat_score_df
