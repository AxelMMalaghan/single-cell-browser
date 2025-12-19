import os

import argparse

os.environ["NUMBA_THREADING_LAYER"] = "tbb"

from sc_browser.ui.dash_app import create_dash_app

from sc_browser.logging_config import configure_logging


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--server.port', dest='port', type=int,

                        default=int(os.getenv("PORT", "8051")))

    parser.add_argument('--host', default="0.0.0.0")

    parser.add_argument('--debug', action='store_true',

                        default=os.getenv("DEBUG", "0") == "1")

    return parser.parse_args()


configure_logging()

app = create_dash_app()

server = app.server

if __name__ == "__main__":
    args = parse_args()

    app.run(host=args.host, port=args.port, debug=args.debug)
