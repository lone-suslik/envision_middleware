# This is a sample Python script.
import os
import re
from collections import OrderedDict

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
import tiledb

app = FastAPI()
index = "studies"

origins = ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

database_uri = os.path.join("/home/suslik/Documents/programming/envision/backend/middle_layer/latest/database")


def get_uri_for_study(study_id: str):
    """

    :param study_id:
    :return: uri: str
    """

    return os.path.join(database_uri, study_id)


def get_uri_for_contrasts(study_id: str):
    """


    :param study_id:
    :return: Path to the contrast array for the specified study
    """

    return os.path.join(database_uri, study_id, "contrasts")


@app.get("/")
async def root():
    return {"message": "The api is responding"}


@app.get("/studies/")
async def studies():
    # TODO: implement checks for path existing

    res = []
    tiledb.ls(database_uri, lambda obj_path, obj_type: res.append(obj_path))
    res = [os.path.basename(os.path.normpath(x)) for x in res]

    return {"message": tiledb.object_type(database_uri), "res": res}


@app.get("/studies/{study_id}/")
async def study_by_id(study_id: str):
    # TODO implement checks for study path existing
    uri = get_uri_for_contrasts(study_id)

    with tiledb.open(uri, 'r') as A:
        a = tiledb.QueryCondition(expression="logFC > 0")
        schema = A.schema
        contrasts = A[:]
        contrasts = dict(zip(contrasts['contrasts'],
                             contrasts['formula']))

    res = {"contrasts": contrasts}

    return {"message": "GET /studies/{study_id}/: the api is responding",
            "res": res,
            "error": None}


@app.get("/studies/{study_id}/contrasts")
async def contrasts_for_study(study_id):
    return {"message": "GET /studies/{study_id}/contrasts: the api is responding"}


@app.get("/studies/{study_id}/contrasts/{contrast_name}/")
async def contrast_summary(study_id, contrast_name):
    return {"message": "GET /studies/{study_id}/contrasts/{contrast_name}/: the api is responding"}


@app.get("/studies/{study_id}/contrasts/{contrast_name}/geneCount/")
async def n_genes_per_contrast(study_id):
    return {"message": "GET /studies/{study_id}/contrasts/{contrast_name}/geneCount/: the api is responding"}


@app.get("/studies/{study_id}/contrasts/pvalues/")
async def study_by_id(study_id):
    return {"message": "GET /studies/{study_id}/contrasts/pvalues/: the api is responding"}
