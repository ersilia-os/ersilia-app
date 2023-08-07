FROM python:3.10.7-bullseye

WORKDIR .

COPY . .

RUN python -m pip install --upgrade pip
RUN python -m pip install ersilia
RUN python -m pip install git+https://github.com/ersilia-os/bentoml-ersilia.git
RUN python -m pip install streamlit==1.23.1

EXPOSE 8501
CMD ["streamlit", "run", "app/app.py"]
