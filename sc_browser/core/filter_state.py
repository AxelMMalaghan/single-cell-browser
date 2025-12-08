from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional

@dataclass
class FilterState:
    """
    Represents the current user selection/filters.

    Fields:

    - genes: List of gene identifiers selected by user.
    - clusters: List of cluster labels selected by user.
    - conditions: List of conditions selected by user.

    - merge_genes: If True, treats multiple genes as a combined signal instead of plotting seperately
    - split_by_condition: If True, views should split their visualization by condition
    - color_scale: Name of color scale used when rendering continuous values

    """

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


@dataclass
class FilterProfile:
    """
    Represents the widget dependencies for different views.

    Fields:

    :param clusters: the cluster widget
    :param : the condition widget
    - conditions: the condition widget
    - samples: the sample widget
    - cell_types: the cell type widget
    - embedding: the embedding widget (for embedding type - projections plots)
    """
    clusters: bool = False
    conditions: bool = False
    samples: bool = False
    cell_types: bool = False
    genes: bool = False
    embedding: bool = False
    split_by_condition: bool = False
    is_3d: bool = False
