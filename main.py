# This is a sample Python script.
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

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


@app.get("/")
async def root():
    return {"message": "The api is responding"}


@app.get("/studies/")
async def studies(prefix="gsf"):
    return {"message": "GET studies: the api is responding"}


@app.get("/studies/{study_id}/")
async def study_by_id(study_id):
    return {"message": "GET /studies/{study_id}/: the api is responding"}


@app.get("/studies/{study_id}/test")
async def test(study_id):
    return {"message": "GET /studies/{study_id}/test: the api is responding"}


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
