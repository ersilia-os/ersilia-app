import streamlit as st
import os
import csv
from ersilia import ErsiliaModel
from ersilia.hub.fetch.fetch import ModelFetcher

ROOT = os.path.dirname(os.path.abspath(__file__))

# Fetch Model
params = st.experimental_get_query_params()
if not "model_id" in params:
    st.error("You need to enter a model identifier as part of the URL, for example: http://localhost:8503/?model_id=eos7yti")
model_id = params["model_id"][0]

mf = ModelFetcher(force_from_hosted=True, hosted_url=None)
if not mf.exists(model_id):
    mf.fetch(model_id)

#get info to populate page
em = ErsiliaModel(model=model_id)
info = em.info()

# Extract the desired values
identifier = info["metadata"]["Identifier"]
slug = info["metadata"]["Slug"]
title = info["metadata"]["Title"]
description = info["metadata"]["Description"]
task = info["metadata"]["Task"]
interpretation = info["metadata"]["Interpretation"]
source_code = info["metadata"]["Source Code"]
publication = info["metadata"]["Publication"]
license = info["metadata"]["License"]

# Theming
st.set_page_config(
    page_title=title,
    layout="wide",
    initial_sidebar_state="expanded"
)
st.markdown('<style>body { color: #50285a; }</style>', unsafe_allow_html=True)


st.title(title,)
st.write(description,)

# Side bar

st.sidebar.title(slug)
st.sidebar.header("Model ID")
st.sidebar.markdown(identifier)
st.sidebar.header("Results")
st.sidebar.markdown(interpretation)
st.sidebar.header("Publication")
st.sidebar.markdown(publication)
st.sidebar.header("Source Code")
st.sidebar.markdown(source_code)
st.sidebar.header("License")
st.sidebar.markdown(license)

# Input
st.header("Input molecules")
st.markdown("Input molecules as a list of SMILES strings. For example:")
smiles = []
with open(os.path.join(ROOT, "..", "data", "example.csv"), "r") as f:
    reader = csv.reader(f)
    for r in reader:
        smiles += [r[0]]
example_smi = ("\n".join(smiles))
st.text(example_smi)

input_molecules = st.text_area(label="",height=125, label_visibility="collapsed")
input_molecules = input_molecules.split("\n")
input_molecules = [inp for inp in input_molecules if inp != ""]

button=st.button('Run')
if (button==True):
  
    def is_valid_input_molecules():
        if len(input_molecules) == 0:
            return False
        for input_molecule in input_molecules:
            mol = Chem.MolFromSmiles(input_molecule)
            if mol is None:
                st.error("Input {0} is not a valid SMILES".format(input_molecule))
                return False
        return True

    if is_valid_input_molecules():

        em.serve()
        df = em.run(input=input_molecules, output="pandas")
        em.close()

        st.dataframe(df, hide_index=True)

        csv_data = df.to_csv(index=False).encode()
        st.download_button(
            "Download as CSV", csv_data, "{}_predictions.csv".format(model_id), "text/csv", key="download-csv"
        )