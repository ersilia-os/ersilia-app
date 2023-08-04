import streamlit as st
import os
import csv
from ersilia import ErsiliaModel
from ersilia.hub.fetch.fetch import ModelFetcher
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import pandas as pd


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

header_color = "#50285a"  # Header background color
separator_color = "#d3d3d3"  # Separator color
header_font_size = "20px"  # Font size for headers

# Theming
st.set_page_config(
    page_title=title,
    layout="wide",
    initial_sidebar_state="expanded"
)
st.markdown(f'<style>body {{ color: {header_color}; background-color: white; }}</style>',unsafe_allow_html=True,)

icon= "🧪"
st.title(f"{icon} {title}")

# Side bar

# Side bar
st.sidebar.title('📊 Model Information')
st.sidebar.title("Description")
st.sidebar.markdown(description)
st.sidebar.header("Identifiers: ")
st.sidebar.markdown(identifier)
st.sidebar.markdown(slug)
st.sidebar.header("Results")
st.sidebar.markdown(interpretation)
st.sidebar.header("Publication")
st.sidebar.markdown(publication)
st.sidebar.header("Source Code")
st.sidebar.markdown(source_code)
st.sidebar.header("License")
st.sidebar.markdown(license)

# Input
st.subheader("Input molecules")
st.markdown("Enter a list of molecules using SMILE notation and each molecule on a separate line, For example:")
smiles = []
with open(os.path.join(ROOT, "..", "data", "example.csv"), "r") as f:
    reader = csv.reader(f)
    for r in reader:
        smiles += [r[0]]
example_smi = ("\n".join(smiles))
st.text(example_smi)

input_molecules = st.text_area(label="",height=125, label_visibility="collapsed")
file_csv= st.file_uploader("Upload csv file", type= ["csv"])
input_molecules = input_molecules.split("\n")
input_molecules = [inp for inp in input_molecules if inp != ""]



button=st.button('Run')
if (button==True):

    if file_csv is not None:
        data_file=pd.read_csv(file_csv)
        input_molecules=data_file['smiles'].tolist()

    def get_molecule_image(smiles):
        m = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(m)
        opts = Draw.DrawingOptions()
        opts.bgColor = None
        im = Draw.MolToImage(m, size=(200, 200), options=opts)
        return im   

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

        st.subheader("Results")

        for idx, r in df.iterrows():
            st.markdown(f"<h4 style='font-size: {header_font_size}; margin-bottom: 10px;'>{idx + 1}: <code>{r['input']}</code></h4>",unsafe_allow_html=True,)
            col1, col2 = st.columns([1, 2])
            image = get_molecule_image(r["input"])
            col1.image(image)

            for col_name, col_value in r.items():
                if col_name != "input" and col_name != "key":
                    text = f"<p style='font-size: 18px; color: white;'>{col_name} : {col_value}</p> "
                    col2.markdown(text, unsafe_allow_html=True)
                
            # Add a separator line with color
            st.markdown(f'<hr style="border: 2px solid {separator_color};">',unsafe_allow_html=True,)

        csv_data = df.to_csv(index=False).encode()
        st.download_button(
            "Download as CSV", csv_data, "{}_predictions.csv".format(model_id), "text/csv", key="download-csv"
        )

