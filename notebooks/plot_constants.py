import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"

tab10 = sns.color_palette("tab10")

# DNA = blue
DNA_COLOR = tab10[0]

# Protein = orange
PROTEIN_COLOR = tab10[1]

# Dayhoff = green
DAYHOFF_COLOR = tab10[2]

PEPTIDE_MOLTYPE_PALETTE = dict(protein=PROTEIN_COLOR, dayhoff=DAYHOFF_COLOR)

PEPTIDE_MOLTYPE_ORDER = "protein", "dayhoff"

PEPTIDE_ALPHABET_KSIZES = {
    "protein": [
        5,
        6,
        7,
        8,
        9,
        10,
    ],
    "dayhoff": [
        15,
        16,
        17,
        18,
        19,
        20,
    ],
}

PEPTIDE_ALPHABET_PALETTES = {
    "protein": sns.dark_palette(PROTEIN_COLOR),
    "dayhoff": sns.dark_palette(DAYHOFF_COLOR),
}
