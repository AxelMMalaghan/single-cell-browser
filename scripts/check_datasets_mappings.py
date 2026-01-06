import json
from pathlib import Path

import anndata as ad

# Adjust these to the absolute or relative paths on your machine
# For your environment: /Users/axelm/PycharmProjects/single-cell-browser/
BASE_DIR = Path(__file__).parent.parent
CONFIG_DIR = BASE_DIR / "config" / "datasets"
DATA_DIR = BASE_DIR / "data"


def check_dataset_mappings():
    print(f"{'DATASET':<30} | {'FIELD':<15} | {'MAPPED TO':<20} | {'STATUS'}")
    print("-" * 85)

    if not CONFIG_DIR.exists():
        print(f"Error: Configuration directory not found at {CONFIG_DIR}")
        return

    config_files = list(CONFIG_DIR.glob("*.json"))
    if not config_files:
        print(f"No .json files found in {CONFIG_DIR}")
        return

    for config_file in config_files:
        with open(config_file, "r") as f:
            try:
                config = json.load(f)
            except json.JSONDecodeError:
                print(f"{config_file.name:<30} | Error: Invalid JSON")
                continue

        dataset_name = config.get("name", config_file.stem)

        # Use the 'path' key from your config files
        rel_h5ad_path = config.get("path", "")
        h5ad_path = (BASE_DIR / rel_h5ad_path).resolve()

        if not h5ad_path.exists():
            print(f"{dataset_name:<30} | Error: .h5ad not found at {h5ad_path}")
            continue

        try:
            # Load in backed mode to avoid memory issues
            adata = ad.read_h5ad(h5ad_path, backed="r")
            obs_cols = adata.obs.columns.tolist()
            obsm_keys = list(adata.obsm.keys())

            # Map the fields as defined in your sc_browser.config.model
            mappings = {
                "cluster": config.get("obs_columns", {}).get("cluster"),
                "condition": config.get("obs_columns", {}).get("condition"),
                "sample": config.get("obs_columns", {}).get("sample"),
                "cell_type": config.get("obs_columns", {}).get("cell_type"),
                "embedding": config.get("embedding_key"),
            }

            for field, mapped_col in mappings.items():
                if not mapped_col:
                    print(
                        f"{dataset_name:<30} | {field:<15} | {'[Not Set]':<20} | SKIP"
                    )
                    continue

                if field == "embedding":
                    exists = mapped_col in obsm_keys
                    container = ".obsm"
                else:
                    exists = mapped_col in obs_cols
                    container = ".obs"

                status = "✅ OK" if exists else f"❌ MISSING in {container}"
                print(f"{dataset_name:<30} | {field:<15} | {mapped_col:<20} | {status}")

            # Specifically for Volcano plots: ensure condition has >1 unique value
            if mappings["condition"] in obs_cols:
                n_unique = adata.obs[mappings["condition"]].nunique()
                if n_unique < 2:
                    print(
                        f"{dataset_name:<30} | {'volcano_check':<15} | {n_unique} groups{'':<11} | ⚠️ NEED >1 FOR DE"
                    )

        except Exception as e:
            print(f"{dataset_name:<30} | Error checking {h5ad_path.name}: {e}")


if __name__ == "__main__":
    check_dataset_mappings()
