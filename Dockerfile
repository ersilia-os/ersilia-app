FROM python:3.10.7-bullseye

WORKDIR .

COPY . .

RUN python -m pip install --upgrade pip
RUN python -m pip install rdkit==2024.3.5
RUN python -m pip install streamlit==1.38.0
RUN python -m pip install ersilia==0.1.36

EXPOSE 8501
CMD ["streamlit", "run", "app/app.py"]
