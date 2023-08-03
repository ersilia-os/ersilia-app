import streamlit as st
from ersilia import ErsiliaModel
from ersilia.hub.fetch.fetch import ModelFetcher

# http://localhost:8503/?model_id=eos7yti

params = st.experimental_get_query_params()
if not "model_id" in params:
    st.error("You need to enter a model identifier as part of the URL, for example: http://localhost:8503/?model_id=eos7yti")
model_id = params["model_id"][0]

mf = ModelFetcher(force_from_hosted=True, hosted_url=None)
if not mf.exists(model_id):
    mf.fetch(model_id)

em = ErsiliaModel(model=model_id)
info = em.info()
st.write(info)
em.serve()
df = em.run(input=["CCCCOCCCC"], output="pandas")
em.close()
st.write(df)

