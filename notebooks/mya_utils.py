from io import StringIO

import pandas as pd


def clean_common_names_of_species(species_name):
    species_name = species_name.replace("_", " ")
    if species_name == "mouse":
        return "house mouse"
    elif species_name == "shrew":
        return "common shrew"
    else:
        return species_name


# http://timetree.org/
# "Estimated time" from human + species

distance_from_human_mya = pd.Series(
    {
        # mus_musculus
        "house mouse": 90,
        # oryctolagus_cuniculus
        "rabbit": 90,
        # rhinolophus_sinicus
        "bat": 96,
        # erinaceus_europaeus
        "hedgehog": 96,
        # sorex_araneus, common shrew
        "common shrew": 96,
        # ornithorhynchus_anatinus
        "platypus": 177,
        # aotus_nancymaae
        "night monkey": 43.2,
        # ceratotherium_simum_simum
        "rhino": 96,
        # camelus_bactrianus
        "camel": 96,
        # lipotes_vexillifer
        "baiji": 96,
        # chinchilla_lanigera
        "chinchilla": 90,
        # macaca_mulatta
        "macaque": 29.4,
        # homo_sapiens
        "human": 0,
        # tupaia_chinensis
        "tupaia": 82,
        # peromyscus_maniculatus_bairdii
        "deer mouse": 90,
        # capra_hircus
        "goat": 96,
        # nannospalax_galili
        "spalax": 90,
        # phascolarctos_cinereus
        "koala": 159,
    }
)
distance_from_human_mya = distance_from_human_mya.sort_values()

# Hard-code in the order from dendrogram
distance_from_human_mya = distance_from_human_mya[
    [
        "human",
        "macaque",
        "night monkey",
        "tupaia",
        "house mouse",
        "deer mouse",
        "spalax",
        "chinchilla",
        "rabbit",
        "hedgehog",
        "common shrew",
        "bat",
        "rhino",
        "camel",
        "baiji",
        "goat",
        "koala",
        "platypus",
    ]
]
distance_from_human_mya


s = """common_name	scientific_name
human	Homo sapiens
house mouse	Mus musculus
spalax	Nannospalax galili
baiji	Lipotes vexillifer
deer mouse	Peromyscus Maniculatus bairdii
tupaia	Tupaia chinensis
bat	Rhinolophus sinicus
rabbit	Oryctolagus cuniculus
hedgehog	Erinaceus europaeus
common shrew	Sorex araneus
platypus	Ornithorhynchus anatinus
night monkey	Aotus nancymaae
rhino	Ceratotherium simum simum
chinchilla	Chinchilla lanigera
koala	Phascolarctos cinereus
macaque	Macaca mulatta
camel	Camelus bactrianus
goat	Capra hircus"""

busco_mammalia_species = pd.read_csv(StringIO(s), sep="\t")
busco_mammalia_species[
    "scientific_lower"
] = busco_mammalia_species.scientific_name.str.lower().str.replace(" ", "_")
busco_mammalia_species = busco_mammalia_species.set_index("common_name")
busco_mammalia_species["mya"] = distance_from_human_mya
busco_mammalia_species = busco_mammalia_species.reset_index()
BUSCO_MAMMALIA_SPECIES = busco_mammalia_species.set_index("scientific_lower")


MYA_ORDER = sorted(distance_from_human_mya.unique())
MYA_COLOR_KWARGS = dict(
        hue='mya',
    hue_order=MYA_ORDER,
    palette='cividis',
)