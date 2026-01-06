import json
from pathlib import Path
from typing import Any, Dict, Optional

from sc_browser.core.dataset import Dataset


def save_dataset_config(
    dataset: Dataset,
    config_root: Path,
    filename: Optional[str] = None,
    schema: str = "anndata_mapped",
) -> Path:
    """
    Persists a Dataset's mapping to a per-dataset JSON file

    Output location: config_root/datasets/<slug>.json
    Schema matches the agreed format

    :param dataset: the new dataset
    :param config_root: the root of the config file
    :param filename: the filename
    :param schema: the schema of the dataset

    """

    datasets_dir = (
        config_root / "datasets"
    )  # TODO:  write to fileshare instead of locally
    datasets_dir.mkdir(parents=True, exist_ok=True)

    if filename is None:
        slug = dataset.name.strip().replace(" ", "_").replace("/", "_")
        filename = f"{slug}.json"

    output_path = datasets_dir / filename

    if dataset.file_path is None:
        raise RuntimeError(
            f"Dataset {dataset.name} is missing a file path, cannot save config"
        )

    project_root = config_root.parent.resolve()
    abs_file = dataset.file_path.resolve()

    try:
        rel_path = abs_file.relative_to(project_root).as_posix()
    except ValueError:
        rel_path = str(abs_file)

    obs_cols: Dict[str, str] = dict(dataset.obs_columns)
    obs_cols.setdefault("cluster", dataset.cluster_key)
    obs_cols.setdefault("condition", dataset.condition_key)

    raw: Dict[str, Any] = {
        "schema": schema,
        "name": dataset.name,
        "group": dataset.group,
        "path": rel_path,
        "obs_columns": obs_cols,
        "embedding_key": dataset.embedding_key,
    }

    with output_path.open("w") as f:
        json.dump(raw, f, indent=2)

    return output_path
