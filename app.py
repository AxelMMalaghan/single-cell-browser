from sc_browser.ui.dash_app import create_dash_app
from sc_browser.logging_config import configure_logging

def main():

    configure_logging()

    app = create_dash_app()
    app.run(debug=True)

if __name__ == '__main__':
    main()