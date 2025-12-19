import os

os.environ["NUMBA_THREADING_LAYER"] = "tbb"

from sc_browser.ui.dash_app import create_dash_app
from sc_browser.logging_config import configure_logging

configure_logging()

app = create_dash_app()
server = app.server

if __name__ == "__main__":
    app.run()