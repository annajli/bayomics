import pandas as pd
import numpy as np

import os
import json
import requests

"""
Encoding module.

Takes the abstract records returned by pubtator_api_module.get_abstract_text()
(list[dict] with "pmid", "title", "abstract") and turns them into dense vector
embeddings for downstream similarity comparison.

MODEL = PubMedBERT

Uses sentence-transformers (local, no API key required) with a biomedical model
by default so domain terminology embeds well. Swap MODEL_NAME for any other
sentence-transformers-compatible model if you'd rather use a general-purpose one.

Install: pip install sentence-transformers >> DONE
"""

MODEL_NAME = "pritamdeka/S-PubMedBert-MS-MARCO"  # biomedical sentence embedding model

_model = None  # lazy-loaded singleton so we don't reload weights on every call


def _get_model():
    global _model
    if _model is None:
        from sentence_transformers import SentenceTransformer

        _model = SentenceTransformer(MODEL_NAME)
    return _model


def build_corpus(abstract_records):
    """
    Combine title + abstract into a single string per record -- this is the
    actual text that gets encoded.

    Parameters
    ----------
    abstract_records : list[dict]
        Output of pubtator_api_module.get_abstract_text().

    Returns
    -------
    list[dict]
        Same records, with an added "text" key (title + abstract concatenated).
        Records with no abstract text are dropped.
    """
    corpus = []
    for record in abstract_records:
        abstract = (record.get("abstract") or "").strip()
        if not abstract:
            continue
        title = (record.get("title") or "").strip()
        combined_text = f"{title}. {abstract}" if title else abstract
        corpus.append({**record, "text": combined_text})
    return corpus


def encode_texts(texts):
    """
    Encode a list of strings into embedding vectors.

    Parameters
    ----------
    texts : list[str]

    Returns
    -------
    np.ndarray
        Shape (len(texts), embedding_dim).
    """
    model = _get_model()
    return model.encode(texts, convert_to_numpy=True, show_progress_bar=False)


def encode_corpus(corpus):
    """
    Encode a corpus built by build_corpus().

    Parameters
    ----------
    corpus : list[dict]
        Each dict must have a "text" key.

    Returns
    -------
    list[dict]
        Same records, with an added "embedding" key (np.ndarray).
    """
    if not corpus:
        return []
    texts = [record["text"] for record in corpus]
    embeddings = encode_texts(texts)
    for record, embedding in zip(corpus, embeddings):
        record["embedding"] = embedding
    return corpus


def cosine_similarity(vec_a, vec_b):
    """Cosine similarity between two 1-D vectors."""
    a = np.asarray(vec_a)
    b = np.asarray(vec_b)
    denom = (np.linalg.norm(a) * np.linalg.norm(b))
    if denom == 0:
        return 0.0
    return float(np.dot(a, b) / denom)


def rank_by_relevance(query_text, encoded_corpus, top_k=5):
    """
    Rank encoded abstracts by cosine similarity to a query string -- e.g. a short
    sentence describing the bnlearn edge ("doxorubicin associated with neoplasms").

    Parameters
    ----------
    query_text : str
    encoded_corpus : list[dict]
        Output of encode_corpus() (each record has "embedding").
    top_k : int
        Number of top-ranked records to return.

    Returns
    -------
    list[dict]
        Top-k records (subset of encoded_corpus), each with an added
        "similarity" score, sorted descending by similarity.
    """
    if not encoded_corpus:
        return []

    query_embedding = encode_texts([query_text])[0]

    scored = []
    for record in encoded_corpus:
        similarity = cosine_similarity(query_embedding, record["embedding"])
        scored.append({**record, "similarity": similarity})

    scored.sort(key=lambda r: r["similarity"], reverse=True)
    return scored[:top_k]