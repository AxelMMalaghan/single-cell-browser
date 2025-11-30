from sc_browser.core.exceptions import ScBrowserError

class DEConfigError(ScBrowserError):
    """DE configuration invalud (bad groupby/group names, etc..)"""
    pass

class DERuntimeError(ScBrowserError):
    """DE failed at runtime (Scanpy error, NaNs, etc)"""
    pass
