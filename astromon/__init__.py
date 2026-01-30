import ska_helpers

__version__ = ska_helpers.get_version(__package__)

from .cross_match import compute_cross_matches
from .db import get_cross_matches

logger = ska_helpers.logging.basic_logger(
    "astromon", level="CRITICAL", format="%(asctime)s %(funcName)-25s: %(message)s"
)


def test(*args, **kwargs):
    """
    Run py.test unit tests.
    """
    import testr

    return testr.test(*args, **kwargs)
