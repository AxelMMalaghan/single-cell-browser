from __future__ import annotations
from dataclasses import dataclass, field
from typing import List

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
    genes: List[str] = field(default_factory=list)
    clusters: List[str] = field(default_factory=list)
    conditions: List[str] = field(default_factory=list)


    merge_genes: bool = False
    split_by_condition: bool = False
    color_scale: str = "viridis"

