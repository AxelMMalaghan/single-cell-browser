from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import Optional

from ..core.state import FilterState

class ExportFormat(Enum):
    PNG = auto()
    SVG = auto()
    PDF = auto()
    JSON = auto()

@dataclass(frozen=True)
class RenderOptions:
    width_pixel: Optional[int] = None
    height_pixel: Optional[int] = None
    dpi: Optional[int] = None
    scale: Optional[float] = None

@dataclass(frozen=True)
class GraphExportRequest:
    """
    Core use case input:

    i.e 'export this view with these filters in this format'
    """
    view_id: str
    filters: FilterState
    format: ExportFormat
    render_options: Optional[RenderOptions] = None
    filename: Optional[str] = None

@dataclass(frozen=True)
class RenderedGraph:
    """
    Output of the rendering backend, still in memory.
    """
    content: bytes
    content_type: str       # e.g. "image/png"
    extension: str         # e.g "png"

@dataclass(frozen=True)
class ExportMetaData:
    """
    Metadata for view, format, filename
    """
    view_id: str
    format: ExportFormat
    filename: str

@dataclass(frozen=True)
class ExportedGraph:
    """
    Exported graph as seen by users
    """
    view_id: str
    format: ExportFormat
    filename: str
    content: Optional[bytes] = None
    uri: Optional[str] = None
    local_path: Optional[Path] = None