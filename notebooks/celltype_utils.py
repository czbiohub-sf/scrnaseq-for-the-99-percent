
import pandas as pd
import scanpy as sc

from path_constants import H5AD

SHARED_CELLTYPES = [
    "Capillary",
    "Alveolar Epithelial Type 2",
    "B cell",
    "T cell",
    "Natural Killer T cell",
    "Macrophage",
    "Monocyte",
    "Dendritic",
    "Fibroblast",
    "Smooth Muscle and Myofibroblast",
]

## Broad to compartment

BROAD_TO_COMPARTMENT = pd.Series({
    "Capillary": "endothelial",
    "Artery": "endothelial",
    "Vein": "endothelial",
    "Alveolar Epithelial Type 2": "epithelial",
    "T cell": "lymphoid",
    "B cell": "lymphoid",
    "Natural Killer T cell": "lymphoid",
    "Natural Killer": "lymphoid",
    "Monocyte": "myeloid",
    "Macrophage": "myeloid",
    "Dendritic": "myeloid",
    "Neutrophil": "myeloid",
    "Fibroblast": "stromal",
    "Smooth Muscle and Myofibroblast": "stromal",
})
broad_to_compartment = BROAD_TO_COMPARTMENT[SHARED_CELLTYPES]
broad_to_compartment


def get_shared_adata():
    adata = sc.read(H5AD)
    adata.obs = adata.obs.reset_index().set_index('cell_id')

    adata_shared = adata[adata.obs.broad_group.isin(SHARED_CELLTYPES)]
    return adata_shared