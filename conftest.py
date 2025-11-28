import logging


def pytest_configure():
    # Global level = WARNING (blocks matplotlib etc.)
    logging.getLogger().setLevel(logging.WARNING)

    # Your package logs at DEBUG
    logging.getLogger("pikachu").setLevel(logging.DEBUG)