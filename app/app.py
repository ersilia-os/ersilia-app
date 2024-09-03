import streamlit as st
import os
import csv
from rdkit import Chem
from ersilia_client import ErsiliaClient
import pandas as pd
import requests

ROOT = os.path.dirname(os.path.abspath(__file__))

st.set_page_config(
    page_title="Ersilia Model Hub App",
    page_icon = os.path.join(ROOT, "..", "data", "Symbol_Plum.png"),
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        'Get Help': None,
        'Report a bug': None,
        'About': "# Ersilia Open Source Initiative. [Read more](https://ersilia.io/about-us) about us. [Support](https://ersilia.io/donate) our mission"
    },
    )


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

def fetch_model_json(url):
    try:
        response = requests.get(url)
        response.raise_for_status()  # Will raise an HTTPError for bad responses
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching the JSON data: {e}")
        return None

def find_model_host_url(json_data, model_id):
    for model in json_data:
        if model.get("Identifier") == model_id:
            return model.get("Host URL")

# Fetch Model URL
try:
    params = st.query_params
    model_id = params["model_id"]
    json_url = "https://ersilia-model-hub.s3.eu-central-1.amazonaws.com/models.json"
    json_data = fetch_model_json(json_url)
    if json_data:
        try:
            host_url = find_model_host_url(json_data, model_id)
        except:
            st.error("Model not found in our database")
    if host_url is None:
        st.error("Model not hosted online")
        
except KeyError as e:
    st.error("You need to enter a model identifier as part of the URL, for example: http://localhost:8500/?model_id=eos7yti")
    exit()


    
client = ErsiliaClient(host_url)
info = client._info()
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
st.sidebar.image(os.path.join(ROOT, "..","data", "Ersilia_Brand.png"), width=150)
st.sidebar.title('Model Information')
st.sidebar.header("Description")
st.sidebar.markdown(description)
st.sidebar.header("Identifiers")
st.sidebar.markdown("{0} | {1}".format(identifier,slug))
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
    written_input = st.text_area(label="input",height=50, label_visibility="collapsed")
    submitted_written = st.form_submit_button("Run")
if (submitted_written==True):
    input_molecules = written_input.split("\n")
    input_molecules = [inp for inp in input_molecules if inp != ""]
    
st.markdown("Or upload a CSV file with a single column named SMILES")
with st.form("csv uploader", clear_on_submit=True):
    file_csv= st.file_uploader(label="input csv", type= ["csv"], label_visibility="collapsed")
    submitted_csv = st.form_submit_button("Run")
    error_placeholder = st.empty()  # Placeholder for error message 
if (submitted_csv==True):
    error_placeholder.empty()
    try:
        data_file=pd.read_csv(file_csv)
        data_file.columns = map(str.upper, data_file.columns)
        if "SMILES" not in data_file.columns:
            error_placeholder.error("Error: The uploaded file must contain a column named 'SMILES'")
        else:
            input_molecules=data_file['SMILES'].tolist()
    except Exception as e:
        error_placeholder.error("An error occurred while processing the file, please upload a .csv file")

batch_size=20
if submitted_written | submitted_csv == True:
    if is_valid_input_molecules():
        with st.spinner('Running the model...'):
            dfs = []
            for n,i in enumerate(range(0, len(input_molecules), batch_size)):
                input_ = input_molecules[i:i+batch_size]
                print(i, len(input_))
                st.toast(f"Calculating molecule batch {n}")
                df = client.run(input_)
                dfs+=[df]
            df_all = pd.concat(dfs)
            st.subheader("Results")
            df_all.rename(columns={"key":"InChiKey", "input": "SMILES"}, inplace=True)
            st.dataframe(df_all, hide_index=True)
            csv_data = df_all.to_csv(index=False).encode()
            st.download_button(
                "Download as CSV", csv_data, "{}_predictions.csv".format(model_id), "text/csv", key="download-csv"
            )


ft = """
<style>
a:link , a:visited{
color: #50285a;  /* theme's text color hex code at 75 percent brightness*/
background-color: transparent;
text-decoration: none;
}

a:hover,  a:active {
color: #BEE6B4; /* theme's primary color*/
background-color: transparent;
text-decoration: underline;
}

#page-container {
  position: relative;
  min-height: 1vh;
}

footer{
    visibility:hidden;
}

.footer {
position: relative;
left: 0;
top:0px;
bottom: 0px;
width: 100%;
background-color: transparent;
color: #50285a; 
text-align: left; 
}
</style>

<div id="page-container">

<div class="footer">
<p style='font-size: 0.9em;'><a style='display: inline; text-align: left;' href="https://ersilia.io/model-hub" target="_blank">Try out other models</a><br 'style= top:3px;'>
<a style='display: inline; text-align: left;' href="https://ersilia.io/" target="_blank">Find more about Ersilia</a><br 'style= top:3px;'>
<a style='display: inline; text-align: left;' href="https://ersilia.io/donate" target="_blank">Support our mission</a></p>
</div>

</div>
"""
st.markdown("---")
st.write(ft, unsafe_allow_html=True)





