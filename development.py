import streamlit as st
from ersilia import ErsiliaModel
from ersilia.hub.fetch.fetch import ModelFetcher

params = st.experimental_get_query_params()
model_id = params["model_id"][0]

mf = ModelFetcher(force_from_hosted=True, hosted_url=None)
if not mf.exists(model_id):
    mf.fetch(model_id)

em = ErsiliaModel(model=model_id)
em.serve()
df = em.run(input=["CCCCOCCCC"], output="pandas")
em.close()
st.write(df)

info = em.info()