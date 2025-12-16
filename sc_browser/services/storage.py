from __future__ import annotations

import shutil
from abc import ABC, abstractmethod
from pathlib import Path
from typing import IO, Iterator, List


class StorageBackend(ABC):
    """
    Abstract interface for file storage (Local, S3, GCS, etc.).
    """

    @abstractmethod
    def write_bytes(self, path: str, data: bytes) -> None:
        pass

    @abstractmethod
    def read_bytes(self, path: str) -> bytes:
        pass

    @abstractmethod
    def list_files(self, prefix: str, suffix: str = "") -> List[str]:
        """List file paths starting with prefix and ending with suffix."""
        pass

    @abstractmethod
    def exists(self, path: str) -> bool:
        pass

    @abstractmethod
    def make_dir(self, path: str) -> None:
        """Ensure a 'directory' exists (no-op on object stores)."""
        pass


class LocalFileSystemStorage(StorageBackend):
    """
    Production-ready local filesystem implementation.
    """

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)

    def _resolve(self, path: str) -> Path:
        # Prevent path traversal attacks
        full_path = (self.root / path).resolve()
        if not str(full_path).startswith(str(self.root)):
            raise ValueError(f"Access denied: {path}")
        return full_path

    def write_bytes(self, path: str, data: bytes) -> None:
        p = self._resolve(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)

    def read_bytes(self, path: str) -> bytes:
        return self._resolve(path).read_bytes()

    def list_files(self, prefix: str, suffix: str = "") -> List[str]:
        # Prefix is treated as a directory in local fs terms
        p = self._resolve(prefix)
        if not p.exists():
            return []

        # Return relative paths as strings
        files = [
            str(f.relative_to(self.root))
            for f in p.glob(f"*{suffix}")
            if f.is_file()
        ]
        return files

    def exists(self, path: str) -> bool:
        return self._resolve(path).exists()

    def make_dir(self, path: str) -> None:
        self._resolve(path).mkdir(parents=True, exist_ok=True)