from __future__ import annotations

from sc_browser.validation.errors import ValidationIssue, ValidationError
from sc_browser.core.dataset import Dataset


def validate_dataset(ds: Dataset) -> None:
    issues: list[ValidationIssue] = []

    obs = ds.adata.obs
    obsm = ds.adata.obsm

    # cluster/condition keys
    if not ds.cluster_key or ds.cluster_key not in obs.columns:
        issues.append(ValidationIssue("DATASET_CLUSTER_KEY", "Missing/invalid cluster_key mapping."))

    if not ds.condition_key or ds.condition_key not in obs.columns:
        issues.append(ValidationIssue("DATASET_CONDITION_KEY", "Missing/invalid condition_key mapping."))

    # embedding key (optional but recommended)
    if ds.embedding_key and ds.embedding_key not in obsm:
        issues.append(
            ValidationIssue(
                "DATASET_EMBEDDING_KEY",
                f"embedding_key '{ds.embedding_key}' not present in adata.obsm.",
            )
        )

    # semantic obs_columns (sample/cell_type optional, but if set must exist)
    sample_key = ds.obs_columns.get("sample")
    if sample_key and sample_key not in obs.columns:
        issues.append(ValidationIssue("DATASET_SAMPLE_KEY", f"sample key '{sample_key}' not in .obs."))

    celltype_key = ds.obs_columns.get("cell_type")
    if celltype_key and celltype_key not in obs.columns:
        issues.append(ValidationIssue("DATASET_CELLTYPE_KEY", f"cell_type key '{celltype_key}' not in .obs."))

    if issues:
        raise ValidationError(issues)
