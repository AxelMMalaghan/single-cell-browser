from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class FilterState:
    """
    Represents the current user selection/filters.

    Fields:

    - genes: List of gene identifiers selected by user.
    - clusters: List of cluster labels selected by user.
    - conditions: List of conditions selected by user.

    - merge_genes: If True, treats multiple genes as a combined signal instead of plotting separately
    - split_by_condition: If True, views should split their visualization by condition
    - color_scale: Name of color scale used when rendering continuous values

    """

    # Global context
    dataset_name: str
    view_id: str

    # Core selections
    genes: List[str] = field(default_factory=list)
    clusters: List[str] = field(default_factory=list)
    conditions: List[str] = field(default_factory=list)

    # New dimensions
    samples: List[str] = field(default_factory=list)
    cell_types: List[str] = field(default_factory=list)

    # Per-request embedding override (e.g. switch PCA/TSNE/UMAP)
    # If None, views should default to Dataset.embedding_key
    embedding: Optional[str] = None

    # Display / plotting options
    split_by_condition: bool = False
    is_3d: bool = False

    # Dimension options
    dim_x: Optional[int] = 0
    dim_y: Optional[int] = 1
    dim_z: Optional[int] = None

    color_scale: Optional[str] = "viridis"

    # To and from dictionary methods

    def to_dict(self) -> Dict[str, Any]:
        return {
            "dataset_name": self.dataset_name,
            "view_id": self.view_id,
            "genes": list(self.genes),
            "clusters": list(self.clusters),
            "conditions": list(self.conditions),
            "samples": list(self.samples),
            "cell_types": list(self.cell_types),
            "embedding": self.embedding,
            "split_by_condition": self.split_by_condition,
            "is_3d": self.is_3d,
            "dim_x": self.dim_x,
            "dim_y": self.dim_y,
            "dim_z": self.dim_z,
            "color_scale": self.color_scale,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> FilterState:
        return cls(
            dataset_name=data.get("dataset_name"),
            view_id=data.get("view_id"),
            genes=list(data.get("genes", [])),
            clusters=list(data.get("clusters", [])),
            conditions=list(data.get("conditions", [])),
            samples=list(data.get("samples", [])),
            cell_types=list(data.get("cell_types", [])),
            embedding=data.get("embedding"),
            split_by_condition=bool(data.get("split_by_condition", False)),
            is_3d=bool(data.get("is_3d", False)),
            dim_x=data.get("dim_x", 0),
            dim_y=data.get("dim_y", 1),
            dim_z=data.get("dim_z", None),
            color_scale=data.get("color_scale", "viridis"),
        )
