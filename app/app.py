import streamlit as st
from rdkit import Chem
import os
import joblib
import pandas as pd
import collections
import csv
import lazyqsar
from rdkit.Chem import AllChem, Draw

# Theming

st.set_page_config(
    page_title="Open Source Malaria Models",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("Open Source Malaria Models")

# Variables

ROOT = os.path.dirname(os.path.abspath(__file__))

model_names = {
    "osm_all_bin1_morgan": "Malaria IC50 (1uM)",
    "osm_all_bin25_morgan": "Malaria IC50 (2.5uM)",
}


# Load models

def load_models():
    models = {}
    models_dir = os.path.join(ROOT, "..", "models")
    for fn in os.listdir(models_dir):
        mn = fn.split(".joblib")[0]
        if mn in model_names.keys():
            display_name = model_names[mn]
            models[display_name] = joblib.load(os.path.join(models_dir, fn))
    return models


models = load_models()
print(models)

# Side bar

st.sidebar.title("Lightweight Models")

st.sidebar.markdown("Lightweight models have been trained using the [Lazy QSAR library](https://github.com/ersilia-os/lazy-qsar).")


st.sidebar.header("Plasmodium falciparum")
texts = model_names.values()
st.sidebar.text("\n".join(texts))


# Input

input_molecules = st.text_area(label="Input molecules")

input_molecules = input_molecules.split("\n")
input_molecules = [inp for inp in input_molecules if inp != ""]

def is_valid_input_molecules():
    if len(input_molecules) == 0:
        return False
    for input_molecule in input_molecules:
        mol = Chem.MolFromSmiles(input_molecule)
        if mol is None:
            st.error("Input {0} is not a valid SMILES".format(input_molecule))
            return False
    return True


def get_molecule_image(smiles):
    m = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(m)
    opts = Draw.DrawingOptions()
    opts.bgColor = None
    im = Draw.MolToImage(m, size=(200, 200), options=opts)
    return im


if is_valid_input_molecules():
    results = collections.OrderedDict()
    results["smiles"] = input_molecules
    print(results)
    for k in models:
        v = models[k]
        results[k] = list(v.predict_proba(input_molecules)[:,1])
    data = pd.DataFrame(results)

    def convert_df(df):
        return df.to_csv(index=False).encode("utf-8")

    csv = convert_df(data)

    st.download_button(
        "Download as CSV", csv, "predictions.csv", "text/csv", key="download-csv"
    )


    #Display predictions
    st.header("Input molecules")
    for v in data.iterrows():
        idx = v[0] +1 #Index from 1
        r = v[1]
        st.markdown("### {0}: `{1}`".format(idx, r["smiles"]))
        cols = st.columns(5)
        image = get_molecule_image(r["smiles"])
        cols[0].image(image)
        texts = [
            "Malaria IC50 (1uM) : {0:.3f}".format(r["Malaria IC50 (1uM)"]),
            "Malaria IC50 (2.5uM): {0:.3f}".format(r["Malaria IC50 (2.5uM)"]),
        ]
        cols[1].text("\n".join(texts))

else:
    st.markdown("Input molecules as a list of SMILES strings. For example:")
    smiles = []
    with open(os.path.join(ROOT, "..", "data", "example.csv"), "r") as f:
        reader = csv.reader(f)
        for r in reader:
            smiles += [r[0]]
    st.text("\n".join(smiles))