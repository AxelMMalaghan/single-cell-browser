from dataclasses import dataclass

@dataclass
class FilterProfile:
    """
    Represents the widget dependencies for different views.

    Fields:

    :param clusters: the cluster widget
    :param conditions: the condition widget
    :param samples: the sample widget
    :param cell_types: the cell-type widget
    :param genes: the gene widget
    :param: embedding: the embedding widget (for embedding type - projections plots)
    :param: split_by_condition:
    :param: is_3d: is 3d
    """
    clusters: bool = False
    conditions: bool = False
    samples: bool = False
    cell_types: bool = False
    genes: bool = False
    embedding: bool = False
    split_by_condition: bool = False
    is_3d: bool = False
