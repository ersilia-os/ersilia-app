import streamlit as st
import os
import csv
from ersilia import ErsiliaModel
from ersilia.hub.fetch.fetch import ModelFetcher
from rdkit import Chem
import pandas as pd


ROOT = os.path.dirname(os.path.abspath(__file__))

# Functions
def is_valid_input_molecules():
    if len(input_molecules) == 0:
        return False
    for input_molecule in input_molecules:
        mol = Chem.MolFromSmiles(input_molecule)
        if mol is None:
            st.error("Input {0} is not a valid SMILES".format(input_molecule))
            return False
    return True


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
    page_icon = os.path.join(ROOT, "..", "data", "Symbol_Plum.png"),
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        'Get Help': None,
        'Report a bug': None,
        'About': "# Ersilia Open Source Initiative. [Read more](https://ersilia.io/about-us) about us. [Support](https://ersilia.io/donate) our mission"
    },
    )

css = r'''
    <style>
        [data-testid="stForm"] {border: 0px}
        [data-testid="stToolbar"] {visibility: hidden !important;}
        footer {visibility: hidden !important;}
    </style>
'''
st.markdown(css, unsafe_allow_html=True)

st.title(title)

# Side bar

# Side bar
st.sidebar.image(os.path.join(ROOT, "..","data", "Ersilia_Brand.png"), width=150)
st.sidebar.title('Model Information')
st.sidebar.header("Description")
st.sidebar.markdown(description)
st.sidebar.header("Identifiers:")
st.sidebar.markdown("* {0}\n* {1}".format(identifier, slug))
st.sidebar.header("Results")
st.sidebar.markdown(interpretation)
st.sidebar.header("References")
st.sidebar.markdown("🔗 [Publication]({})".format(publication))
st.sidebar.markdown("🔗 [Source Code]({})".format(source_code))
st.sidebar.header("License")
st.sidebar.markdown(license)


# Input
st.subheader("Input molecules")

st.markdown("Enter a list of molecules using SMILES notation and each molecule on a separate line")
smiles = []
with open(os.path.join(ROOT, "..", "data", "example.csv"), "r") as f:
    reader = csv.reader(f)
    for r in reader:
        smiles += [r[0]]
example_smi = ("\n".join(smiles))
st.text(example_smi)
with st.form("text uploader", clear_on_submit=True):
    written_input = st.text_area(label="",height=50, label_visibility="collapsed")
    submitted_written = st.form_submit_button("Run")
if (submitted_written==True):
    input_molecules = written_input.split("\n")
    input_molecules = [inp for inp in input_molecules if inp != ""]
    
st.markdown("Or upload a CSV file with a single column named SMILES")
with st.form("csv uploader", clear_on_submit=True):
    file_csv= st.file_uploader(label="", type= ["csv"], label_visibility="collapsed")
    submitted_csv = st.form_submit_button("Run")
if (submitted_csv==True): 
    data_file=pd.read_csv(file_csv)
    input_molecules=data_file['smiles'].tolist()
    
if submitted_written | submitted_csv == True:
    if is_valid_input_molecules():
        with st.spinner('Running the model...'):
            em.serve()
            df = em.run(input=input_molecules, output="pandas")
            em.close()
            st.subheader("Results")
            df.rename(columns={"key":"InChiKey", "input": "SMILES"}, inplace=True)
            st.dataframe(df, hide_index=True)
            csv_data = df.to_csv(index=False).encode()
            st.download_button(
                "Download as CSV", csv_data, "{}_predictions.csv".format(model_id), "text/csv", key="download-csv"
            )