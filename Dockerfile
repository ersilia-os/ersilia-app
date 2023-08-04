FROM python:3.8.11-bullseye

WORKDIR .

COPY . .

RUN python -m pip install --upgrade pip
RUN python -m pip install streamlit
RUN python -m pip install git+https://github.com/ersilia-os/bentoml-ersilia
RUN python -m pip install ersilia

EXPOSE 8501
CMD ["streamlit", "run", "app/app.py"]