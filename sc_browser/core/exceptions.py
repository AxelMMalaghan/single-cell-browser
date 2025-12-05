

class ScBrowserError(Exception):
    """Base exception for all sc_browser errors"""
    pass

class ConfigError(ScBrowserError):
    """Invalid or inconsistent demo1-dataset.json or global config"""
    pass

class DatasetSchemaError(ScBrowserError):
    """
    AnnData / schema doesn't match what Dataset/Adapter expects
    missing obs/obsm keys, wrong dtypes, etc
    """
    pass

class AdapterSelectionError(ScBrowserError):
    """No suitable DatasetAdapter found for a given adata/config"""
    pass
