import ska_helpers

__version__ = ska_helpers.get_version(__package__)

from .db import get_cross_matches  # noqa
from .cross_match import compute_cross_matches  # noqa
