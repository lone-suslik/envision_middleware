#!/usr/bin/env python

from elasticsearch import Elasticsearch
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional, Dict
from loguru import logger
import json
import numpy as np
import orjson
import os
import pandas
import requests
import tiledb

description = """
This is a REST API application that serves as a middleware for 
the envision app.

It provides access to the following databases:
1 - Tiledb database that stores TPM-level expression data and summary statistics
2 - Elasticsearch instance that stores various gene annotation information

The following species are currently supported:

- Canis lupus familiaris
- Homo sapiens
- Gallus gallus
- Mus musculus
- Rattus norvegicus
- Caenorhabditis elegans
- Drosophila melanogaster
- Saccharomyces cerevisiae

    
"""

"""
IMPORTANT: READ BEFORE CHANGING CODE
Coding conventions for this tool:

Tiledb and elasticsearch accessor functions should abstract interactions with respective
databases. 
ALL functions should be located in the dedicated section (indicated by the comments)

NAMING CONVENTIONS:
Functions that interact with tiledb should start with tdb_. 
Functions that interact with elasticsearch should start with es_.
Functions that implement API endpoints should start with api_.
Helper functions that do not fall under the above should start with helper_.

API CONVENTIONS:
GET API functions should start with api_get

POST API functions should start with api_post

API endpoints should be async unless there is a reason for the opposite.

Query parameters in API endpoints are not allowed - use POST with a json payload (this is for security
reasons and consistency with endpoints that require large bodies of data, for example expression slicing endpoints

API endpoints should have 2 decorators - one without trailing slash and one with trailing slash and , include_in_schema=False
This is explained at https://github.com/tiangolo/fastapi/issues/2060

ADDITIONALLY:
Pydantic models should be located in the pydantic model section (indicated by the comments)

DATABASE DESCRIPTION:
The following information is available in the tiledb database:
TODO - describe tiledb schema 
Below is the description of the elasticsearch schema :

reactome pathways:
    "reactome_{species_name}" : {
        ["genes"] : {"type" : "keyword"},
        "pathway_id" : { "type" : "keyword"},
        "pathway_name" : {"type" : "keyword"},
        "pathway_url" : {"type" : "keyword"}
    }

"""

# Defaults and globals
app = FastAPI(
    title="Envision API",
    description=description,
    version="0.0.0.0.0.1*e-16",
    contact={"name": "Petr Volkov"}
)

origins = ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

database_uri = os.path.join("/home/suslik/Documents/programming/envision/backend/middle_layer/latest/database")
es_uri = "http://localhost:9200"


# pydantic models
class StatsQuery(BaseModel):
    """
    Pydantic class to hold conditional stats filter requests
    """
    filter: str


class PathwayQuery(BaseModel):
    """
    Pydantic class to hold queries to search for pathways
    query - a wildcard string for full-text pathway search
    size - number of items to return
    """
    query: str
    size: Optional[int] = 10
    fields: Optional[List[str]] = None


class ExpressionSlice(BaseModel):
    genes: Optional[List[str]] = None
    samples: Optional[List[str]] = None


class GeneIDQuery(BaseModel):
    genes: List[str]


# ES accessor functions
def es_index_search_uri(index: str, filter_path: str = ""):
    res = f"{es_uri}/{index}/_search"

    if filter_path:
        res = f"{res}?filter_path={filter_path}"

    print(res)
    return res


def es_search_request(index: str, payload: Dict, filter_path: str = ""):
    uri = es_index_search_uri(index, filter_path)
    return requests.post(uri, json=payload)


def es_reactome_search_request(index: str, payload: Dict, filter_path: str = ""):
    return es_search_request(f"reactome_{index}", payload, filter_path)


def es_get_all_reactome_pathways(species: str):
    list_pathways_query = {
        "aggs": {
            "unique_pathways": {
                "terms": {
                    "field": "pathway_name", "size": 10000
                }
            }
        },
        "size": 0
    }
    response = es_reactome_search_request(index=species,
                                          payload=list_pathways_query).json()

    response = [x["key"] for x in response["aggregations"]["unique_pathways"]["buckets"]]

    return response


def es_filter_reactome_pathways_by_query(species: str, query: PathwayQuery):
    """
    :param species:
    :param query:
    :return:

    For response filtering (filter_path parameter) see:
    https://www.elastic.co/guide/en/elasticsearch/reference/current/common-options.html#common-options-response-filtering
    """
    payload = {
        "size": query.size,
        "query": {
            "wildcard": {
                "pathway_name": {
                    "value": query.query
                }
            }
        }
    }

    if query.fields:
        payload["fields"] = query.fields
        payload["_source"] = False

    response = es_reactome_search_request(index=species,
                                          payload=payload,
                                          filter_path="hits.hits.fields").json()

    return response


def es_get_ensembl_gene_symbols_by_id(species: str, query: GeneIDQuery):
    payload = {
        "query": {
            "terms": {
                "gene_id": query.genes
            }
        },
        "size": len(query.genes)
    }
    index_name = f"ensembl_gene_{species}"
    response = es_search_request(index_name,
                                 payload=payload,
                                 filter_path="hits.hits").json()["hits"]["hits"]
    response = {r["_source"]["gene_id"]: r["_source"] for r in response}

    return response


# TILEDB accessor functions
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


def tdb_uri_for_tpm_sparse(study_id: str):
    return os.path.join(tdb_uri_for_study(study_id), "tpm")


def tdb_uri_for_tpm_dense(study_id: str):
    return os.path.join(tdb_uri_for_study(study_id), "tpm_dense")


async def tdb_get_contrasts_for_study(study_id: str):
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


async def tdb_get_genes_for_stats_query(study_id: str, contrast_id: str, query: StatsQuery):
    await tdb_verify_contrast_uri(study_id, contrast_id)
    contrast_id = bytes(contrast_id, 'utf-8')
    uri = tdb_uri_for_stats(study_id)

    try:
        with tiledb.open(uri, 'r') as A:
            qc = tiledb.QueryCondition(query.filter)
            response = A.query(attr_cond=qc, use_arrow=False).multi_index[:, contrast_id]
            res = {
                "genes": list(response['genes'])
            }

    except tiledb.TileDBError as err:
        raise HTTPException(status_code=404,
                            detail=err.message)

    return res


# helpers
def helper_flatten_es_response_dict(response: Dict):
    """
    Elastic search returns values for fields as lists,
    even if there is only one value.

    This function transforms something like:
    {
      "pathway_name" : [
         "p53-Independent G1/S DNA damage checkpoint"
      ],
      "pathway_url" : [
         "https://reactome.org/PathwayBrowser/#/R-HSA-69613"
      ]
    }
    into
    {
      "pathway_name" : "p53-Independent G1/S DNA damage checkpoint",
      "pathway_url" : "https://reactome.org/PathwayBrowser/#/R-HSA-69613"
    }

    It will throw ValueError exception if there are fields with more than 1 value.
    Solution borrowed from: https://stackoverflow.com/questions/10756427/loop-through-all-nested-dictionary-values
    :param response:
    :return:
    """
    for key, value in response.items():
        if type(value) is dict:
            response[key] = helper_flatten_es_response_dict(value)
        elif type(value) is list:
            if len(value) != 1:
                raise ValueError
            else:
                response[key] = response[key][0]

    return response


# API
@app.get("/")
async def api_root():
    return {"message": "The api is responding"}


@app.get("/studies")
@app.get("/studies/", include_in_schema=False)
async def api_get_studies():
    # TODO: implement checks for path existing

    res = await tdb_get_studies()

    return {"message": tiledb.object_type(database_uri), "res": res}


@app.get("/studies/{study_id}")
@app.get("/studies/{study_id}/", include_in_schema=False)
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
@app.get("/studies/{study_id}/contrasts/", include_in_schema=False)
async def api_get_contrasts_for_study(study_id: str):
    res = await api_get_study_by_id(study_id)
    return res


@app.get("/studies/{study_id}/contrasts/{contrast_id}")
@app.get("/studies/{study_id}/contrasts/{contrast_id}/", include_in_schema=False)
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


@app.post("/studies/{study_id}/contrasts/{contrast_id}/query")
@app.post("/studies/{study_id}/contrasts/{contrast_id}/query/", include_in_schema=False)
async def api_post_contrast_filter(study_id: str, contrast_id: str, query: StatsQuery):
    """
    :param study_id:
    :param contrast_id:
    :param query:
         Conditional filter to subset genes based on available statistics
    :return:
    """
    return await tdb_get_genes_for_stats_query(study_id, contrast_id, query)


@app.get("/studies/{study_id}/{contrast_id}/query")
@app.get("/studies/{study_id}/{contrast_id}/query/", include_in_schema=False)
async def api_get_study_tpm_by_query(study_id: str, contrast_id: str):
    """
    This is a testing endpoint to try to use queries to filter the datasets
    :param study_id:
    :param contrast_id:
    :return:
    """
    query = {"filter": "pvalue < 0.0001"}
    query = StatsQuery(**query)

    gene_ids = await tdb_get_genes_for_stats_query(study_id, contrast_id, query)
    gene_ids = gene_ids['genes']
    sample_ids = []

    # tdb_splice_tpm(study_id, contrast_id, gene_ids, sample_ids)
    uri = tdb_uri_for_tpm_sparse(study_id)
    with tiledb.open(uri, 'r') as A:
        if sample_ids:
            res = A.df[gene_ids, :]
        else:
            res = A.df[gene_ids, sample_ids]

    res = pandas.pivot(res, index='genes', columns='sample', values='expr')

    return json.loads(res.to_json())


@app.post("/studies/{study_id}/tpm/slice")
@app.post("/studies/{study_id}/tpm/slice/", include_in_schema=False)
async def api_splice_tpm(study_id: str, slice: ExpressionSlice):
    """
    This is the central endpoint to request expression data for a particular table
    :param study_id:
    :param slice:
    :return:
    """
    # TODO refactor
    uri = tdb_uri_for_tpm_sparse(study_id)
    sample_ids = slice.samples
    gene_ids = slice.genes

    with tiledb.open(uri, 'r') as A:
        if gene_ids and sample_ids:
            res = A.df[gene_ids, sample_ids].pivot(index="sample", columns="genes", values="expr")
        elif gene_ids:
            res = A.df[gene_ids, :]
        elif sample_ids:
            res = A.df[:, sample_ids]

    if not (gene_ids or sample_ids):
        uri = tdb_uri_for_tpm_dense(study_id)
        with tiledb.open(uri, 'r') as A:
            res = A.query(attrs=['expr'])[:]  # ['expr']

        res = orjson.dumps(res,
                           option=orjson.OPT_SERIALIZE_NUMPY)

    return res


@app.get("/reactome/{species}/pathways")
@app.get("/reactome/{species}/pathways/", include_in_schema=False)
async def api_get_all_reactome_pathways(species: str):
    return es_get_all_reactome_pathways(species)


@app.post("/reactome/{species}/pathways/query")
@app.post("/reactome/{species}/pathways/query/", include_in_schema=False)
async def api_post_reactome_studies_by_query(species: str, query: PathwayQuery):
    """
    :param species:
    :param query:
    :return:
    """
    response = es_filter_reactome_pathways_by_query(species, query)["hits"]["hits"]
    [helper_flatten_es_response_dict(x) for x in response]

    return response


@app.post("/ensembl/{species}")
@app.post("/ensembl/{species}/", include_in_schema=False)
async def api_get_ensembl_gene_symbols_by_id(species: str, query: GeneIDQuery):
    """

    :param species:
    :param genes:
    :return:
    """
    response = es_get_ensembl_gene_symbols_by_id(species, query)
    return response
