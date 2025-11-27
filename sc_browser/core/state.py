from __future__ import annotations
from dataclasses import dataclass, field
from typing import List

@dataclass
class FilterState:
    """
    Represents the current user selection/filters
    """
    genes: List[str] = field(default_factory=list)
    clusters: List[str] = field(default_factory=list)
    conditions: List[str] = field(default_factory=list)


    merge_genes: bool = False
    split_by_condition: bool = False
    color_scale: str = "viridis"

