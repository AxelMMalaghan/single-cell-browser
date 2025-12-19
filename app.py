import os
import socket

# Set Numba threading before other imports to prevent crashes
os.environ["NUMBA_THREADING_LAYER"] = "tbb"

from sc_browser.ui.dash_app import create_dash_app
from sc_browser.logging_config import configure_logging

configure_logging()

app = create_dash_app()
server = app.server


def find_free_port(start_port: int) -> int:
    """Finds an available port starting from start_port."""
    port = start_port
    while port < start_port + 100:  # Try up to 100 ports
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('localhost', port)) != 0:
                return port
        port += 1
    return start_port


if __name__ == "__main__":
    # 1. Get preferred port from env or default to 8051
    preferred_port = int(os.getenv("PORT", "8051"))

    # 2. Check if the port is free; if not, find the next available one
    final_port = find_free_port(preferred_port)

    debug = os.getenv("DEBUG", "0") == "1"

    if final_port != preferred_port:
        print(f"Warning: Port {preferred_port} was taken. Starting on {final_port}")

    app.run(host="0.0.0.0", port=final_port, debug=debug)