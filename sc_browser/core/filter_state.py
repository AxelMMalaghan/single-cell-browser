from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import List, Optional, Dict, Any


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
    merge_genes: bool = False
    split_by_condition: bool = False
    is_3d: bool = False

    # Dimension options
    dim_x: int = 0
    dim_y: int = 1
    dim_z: int = 2

    color_scale: str = "viridis"

    @classmethod
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

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
            merge_genes=bool(data.get("merge_genes", False)),
            split_by_condition=bool(data.get("split_by_condition", False)),
            is_3d=bool(data.get("is_3d", False)),
            dim_x=data.get("dim_x", 0),
            dim_y=data.get("dim_y", 1),
            dim_z=data.get("dim_z", 2),
            color_scale=data.get("color_scale", "viridis"),
        )


