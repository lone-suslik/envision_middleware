# This is a sample Python script.
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
import numpy as np
import os
import tiledb

app = FastAPI()

origins = ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

database_uri = os.path.join("/home/suslik/Documents/programming/envision/backend/middle_layer/latest/database")


def tdb_uri_for_study(study_id: str):
    """

    :param study_id:
    :return: uri: str
    """

    return os.path.join(database_uri, study_id)


def tdb_uri_for_contrasts(study_id: str):
    """


    :param study_id:
    :return: Path to the contrast array for the specified study
    """

    return os.path.join(tdb_uri_for_study(study_id), "contrasts")


def tdb_uri_for_stats(study_id: str):
    return os.path.join(tdb_uri_for_study(study_id), "stats")


async def tdb_get_contrasts_for_study(study_id):
    uri = tdb_uri_for_contrasts(study_id)
    with tiledb.open(uri, 'r') as A:
        # a = tiledb.QueryCondition(expression="logFC > 0")
        contrasts = A[:]
        contrasts = dict(zip(contrasts['contrasts'],
                             contrasts['formula']))
    return contrasts


async def tdb_get_studies(to_bytes: bool = False):
    res = []
    tiledb.ls(database_uri, lambda obj_path, obj_type: res.append(obj_path))
    res = [os.path.basename(os.path.normpath(x)) for x in res]
    if to_bytes:
        res = [bytes(str(x), 'utf-8') for x in res]

    return res


@app.get("/")
async def root():
    return {"message": "The api is responding"}


@app.get("/studies/")
async def api_get_studies():
    # TODO: implement checks for path existing

    res = await tdb_get_studies()

    return {"message": tiledb.object_type(database_uri), "res": res}


@app.get("/studies/{study_id}/")
async def api_get_study_by_id(study_id: str):
    # check that study can be found
    studies = await tdb_get_studies()
    if study_id not in studies:
        raise HTTPException(status_code=404,
                            detail=f"Study {study_id} not available in the database")
    contrasts = await tdb_get_contrasts_for_study(study_id)
    res = {"contrasts": contrasts}

    return {"message": "GET /studies/{study_id}/: the api is responding",
            "res": res,
            "error": None}


@app.get("/studies/{study_id}/contrasts")
async def api_get_contrasts_for_study(study_id: str):
    res = await api_get_study_by_id(study_id)
    return res


@app.get("/studies/{study_id}/contrasts/{contrast_id}/")
async def api_get_contrast_summary(study_id: str, contrast_id: str):
    await tdb_verify_contrast_uri(study_id, contrast_id)
    contrast_id = bytes(contrast_id, 'utf-8')

    uri = tdb_uri_for_stats(study_id)
    with tiledb.open(uri, 'r') as A:
        pvalues = A[:, contrast_id]['pvalue']
        fdr = A[:, contrast_id]['fdr']

    print(np.where(pvalues < 0.05))

    res = {
        "n_genes": len(pvalues),
        "n_fdr_significant_genes": len(np.where(fdr < 0.05)[0]),
        "n_p_significant_genes": len(np.where(pvalues < 0.05)[0])
    }

    return {"message": "GET /studies/{study_id}/contrasts/{contrast_name}/: the api is responding", "res": res}


@app.post("/studies/{study_id}/contrasts/{contrast_id}/filter")
async def api_post_contrast_filter(study_id: str, contrast_id: str):
    await tdb_verify_contrast_uri(study_id, contrast_id)
    contrast_id = bytes(contrast_id, 'utf-8')


async def tdb_verify_contrast_uri(study_id: str, contrast_id: str):
    # check that study and contrast can be found
    contrast_id = bytes(contrast_id, 'utf-8')
    studies = await tdb_get_studies()
    if study_id not in studies:
        raise HTTPException(status_code=404,
                            detail=f"Study {study_id} not available in the database")
    contrasts = dict(await tdb_get_contrasts_for_study(study_id))
    if contrast_id not in contrasts:
        raise HTTPException(status_code=404,
                            detail=f"Contrast {contrast_id} not available for the study {study_id}")