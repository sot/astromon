import hashlib
import inspect

# import logging
# import functools
import json
from pathlib import Path

from astropy.table import Table


def is_writable(dirpath):
    import os
    import sys
    import uuid

    if sys.platform.startswith("linux"):
        return os.access(dirpath, os.W_OK)

    dirpath = Path(dirpath)
    test_file = dirpath / f"{uuid.uuid4()}.txt"
    try:
        with open(test_file, "w") as test:
            test.write("hello")
            return True
    except OSError:
        return False
    finally:
        if test_file.exists():
            test_file.unlink()


class Cache:
    """
    Decorator to cache function return values.
    """

    def __init__(self, name, fmt="json", dir_attribute=None, cache_location=None):
        self.name = name
        self.dir_attribute = dir_attribute
        self.fmt = fmt
        self.cache_location = Path() if cache_location is None else cache_location

    def __call__(self, func):
        class Wrapper:
            def __init__(self, name, dir_attribute, fmt, cache_location):
                self.name = name
                self.dir_attribute = dir_attribute
                self.fmt = fmt
                self.cache_location = cache_location

            # this implements the decorator for functions
            def __call__(self, *args, **kwargs):
                hsh = argument_hash(func, *args, **kwargs)
                filename = self.cache_location / f"{self.name}::{hsh}"
                return Wrapper.do(filename, self.fmt, *args, **kwargs)

            # this implements clear_cache for functions
            def clear_cache(self, *args, **kwargs):
                hsh = argument_hash(func, *args, **kwargs)
                filename = self.cache_location / f"{self.name}::{hsh}"
                return self.clear_cache_(filename, self.fmt, *args, **kwargs)

            # this implements clear_all for functions
            def clear_all(self):
                for filename in self.cache_location.glob(f"{self.name}::*"):
                    filename.unlink(missing_ok=True)

            def files(self):
                return list(self.cache_location.glob(f"{self.name}::*"))

            def __get__(self, obj, objtype_):
                class Enhanced:
                    # this implements the decorator for class methods
                    def __call__(enh, *args, **kwargs):
                        hsh = argument_hash(func, obj, *args, **kwargs)
                        if self.dir_attribute is None:
                            filename = f"{self.name}::{hsh}"
                        else:
                            filename = (
                                getattr(obj, self.dir_attribute) / f"{self.name}::{hsh}"
                            )
                        return Wrapper.do(filename, self.fmt, obj, *args, **kwargs)

                    # this implements clear_cache for class methods
                    def clear_cache(enh, *args, **kwargs):
                        hsh = argument_hash(func, obj, *args, **kwargs)
                        if self.dir_attribute is None:
                            filename = f"{self.name}::{hsh}"
                        else:
                            filename = (
                                getattr(obj, self.dir_attribute) / f"{self.name}::{hsh}"
                            )
                        return self.clear_cache_(filename, self.fmt, *args, **kwargs)

                    def clear_all(enh):
                        if self.dir_attribute is None:
                            cache_location = self.cache_location
                        else:
                            cache_location = self.cache_location / getattr(
                                obj, self.dir_attribute
                            )
                        for filename in cache_location.glob(f"{self.name}::*"):
                            filename.unlink(missing_ok=True)

                    def files(enh):
                        if self.dir_attribute is None:
                            cache_location = self.cache_location
                        else:
                            cache_location = self.cache_location / getattr(
                                obj, self.dir_attribute
                            )
                            return list(cache_location.glob(f"{self.name}::*"))

                return Enhanced()

            @staticmethod
            def do(filename, fmt, *args, **kwargs):
                filename = FORMATS[fmt].sanitize_filename(filename)
                if filename.exists():
                    return FORMATS[fmt].load(filename)
                result = func(*args, **kwargs)
                FORMATS[fmt].save(result, filename)
                return result

            @staticmethod
            def clear_cache_(filename, fmt, *args, **kwargs):
                filename = FORMATS[fmt].sanitize_filename(filename)
                filename.unlink(missing_ok=True)

        return Wrapper(self.name, self.dir_attribute, self.fmt, self.cache_location)


def argument_hash(func, *args, **kwargs):
    # kwargs are sorted to be consistent. All other arguments are sorted in the signature already.
    kwargs = {k: kwargs[k] for k in sorted(kwargs)}

    sig = inspect.signature(func)
    sig_args = sig.bind(*args, **kwargs)
    sig_args.apply_defaults()
    sig_args = sig_args.arguments.copy()
    sig_args.pop("self", None)

    # this cache treats lists and tuple arguments the same way, so you don't get two entries,
    # one for [1, 2, 3] and one for (1, 2, 3)
    for k, v in sig_args.items():
        if isinstance(v, list):
            sig_args[k] = tuple(v)
    hsh = hashlib.md5(json.dumps(sig_args).encode()).hexdigest()
    return hsh


class TableCache:
    @staticmethod
    def save(table, filename, force=False):
        if filename.exists() and not force:
            return
        if not filename.parent.exists():
            filename.parent.mkdir(parents=True, exist_ok=True)
        if not isinstance(table, Table):
            table = Table(table)
        if "description" in table.meta:
            table.meta["descript"] = table.meta["description"]
            del table.meta["description"]
        table.write(filename, overwrite=force)

    @staticmethod
    def load(filename):
        if filename.exists():
            table = Table.read(filename)
            table.convert_bytestring_to_unicode()
            return table

    @staticmethod
    def sanitize_filename(filename):
        filename = str(filename)
        if not filename.endswith((".fits", ".fits.gz")):
            filename = Path(str(filename) + ".fits.gz")
        return Path(filename)


class JsonCache:
    @staticmethod
    def save(result, filename, force=False):
        if filename.exists() and not force:
            return
        if not filename.parent.exists():
            filename.parent.mkdir(parents=True, exist_ok=True)
        with open(filename, "w") as fh:
            json.dump(result, fh)

    @staticmethod
    def load(filename):
        if filename.exists():
            with open(filename, "r") as fh:
                return json.load(fh)

    @staticmethod
    def sanitize_filename(filename):
        filename = str(filename)
        if not (filename.endswith(".json")):
            filename = filename + ".json"
        return Path(filename)


FORMATS = {
    "table": TableCache,
    "json": JsonCache,
}
