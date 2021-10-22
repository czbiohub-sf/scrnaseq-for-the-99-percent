import matplotlib.pyplot as plt
import seaborn as sns
from umap import UMAP


def plot_umap(transformed, *args, **kwargs):
    x = transformed[:, 0]
    y = transformed[:, 1]

    if 'ax' not in kwargs:
        ax = plt.gca()
    else: 
        ax = kwargs.pop('ax')
    
    # From scanpy
    size = 120000 / transformed.shape[0]
    sns.scatterplot(x, y, s=size, edgecolor='none', ax=ax, *args, **kwargs)
    
    sns.despine(left=True, bottom=True)
    ax.set(xlabel='UMAP 1', ylabel='UMAP 2', xticks=[], yticks=[])

def do_umap(similarities, random_state=0, y=None, *args, **kwargs):
    umapper = UMAP(metric='precomputed', random_state=0, *args, **kwargs)
    
    # Convert similarity to dissimilarity which is approximately a distance
    transformed = umapper.fit_transform(1-similarities, y=y)
    return transformed, umapper

def umap_and_plot(similarities, random_state=0, hue=None, *args, **kwargs):
    transformed, umapper = do_umap(similarities)
    plot_umap(transformed, hue=hue, *args, **kwargs)
    return transformed, umapper
