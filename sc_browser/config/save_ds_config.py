import json
import re
from pathlib import Path
from typing import Any, Dict, Optional


def _slugify(name: str) -> str:
    s = name.strip()
    s = re.sub(r"\s+", "_", s)
    s = s.replace("/", "_")
    s = re.sub(r"[^A-Za-z0-9_.-]+", "", s)
    return s or "dataset"


def write_dataset_config_for_uploaded_file(
    *,
    config_root: Path,
    dataset_name: str,
    file_path: Path,
    group: str = "Imported",
    schema: str = "anndata_mapped",
    embedding_key: Optional[str] = None,
    obs_columns: Optional[Dict[str, str]] = None,
    filename: Optional[str] = None,
) -> Path:
    """
    Create config_root/datasets/<slug>.json *without* needing a Dataset object.
    Uses the SAME key names as save_dataset_config(): schema, name, group, path, obs_columns, embedding_key.
    """
    datasets_dir = config_root / "datasets"
    datasets_dir.mkdir(parents=True, exist_ok=True)

    if filename is None:
        filename = f"{_slugify(dataset_name)}.json"

    output_path = datasets_dir / filename

    project_root = config_root.parent.resolve()
    abs_file = file_path.resolve()

    try:
        rel_path = abs_file.relative_to(project_root).as_posix()
    except ValueError:
        rel_path = str(abs_file)

    cols: Dict[str, str] = dict(obs_columns or {})
    # Don’t set cluster/condition yet; user will do mapping
    # (you can set defaults here if you want)

    raw: Dict[str, Any] = {
        "schema": schema,
        "name": dataset_name,
        "group": group,
        "path": rel_path,
        "obs_columns": cols,
        "embedding_key": embedding_key,
    }

    with output_path.open("w", encoding="utf-8") as f:
        json.dump(raw, f, indent=2)

    return output_path
