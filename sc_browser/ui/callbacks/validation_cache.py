from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class ValidSets:
    clusters: set[str]
    conditions: set[str]
    samples: set[str]
    cell_types: set[str]
    genes: set[str]

def build_valid_sets(ds) -> ValidSets:
    clusters = set(ds.clusters.astype(str).unique())
    conditions = set(ds.conditions.astype(str).unique()) if ds.conditions is not None else set()
    samples = set(ds.samples.astype(str).unique()) if getattr(ds, "samples", None) is not None else set()
    cell_types = set(ds.cell_types.astype(str).unique()) if getattr(ds, "cell_types", None) is not None else set()

    genes = (
        set(getattr(ds, "genes", []))
        or set(getattr(ds, "gene_names", []))
        or set(ds.adata.var_names)
    )

    return ValidSets(
        clusters=clusters,
        conditions=conditions,
        samples=samples,
        cell_types=cell_types,
        genes=genes,
    )
