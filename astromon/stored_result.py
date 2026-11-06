import hashlib
import inspect

# import logging
# import functools
import json
import os
import pickle
import shutil
import sys
import tempfile
import uuid
from pathlib import Path
from typing import Any, Callable, Protocol

from astropy.table import Table


class StorageProtocol(Protocol):
    def glob(self, pattern: str) -> list[Path]: ...
    def path(self, name: str) -> Path: ...
    def work_path(self, name: str) -> Path: ...


class StorableProtocol(Protocol):
    storage: StorageProtocol


class Storage:
    """
    Manages storage of files in a working directory and an optional archive directory.

    The Storage class handles files in two locations: a working directory and an optional archive
    directory. The working directory is where files are created and modified during processing,
    while the archive directory is a long-term storage location. When files are archived, they are
    copied from the working directory to the archive directory.

    The Storage class provides methods to access files in both locations, search for files
    matching a pattern, and archive files based on specified patterns.

    The Storage class DOES NOT take ownership of the directories. It only manages file paths.

    Parameters
    ----------
    workdir : str or Path, optional
        Path to the working directory. If None, the current directory is used.
    archive_dir : str or Path, optional
        Path to the archive directory. If None, archiving is disabled.
    instance : StorableProtocol, optional
        An instance of a class implementing the StorableProtocol. If provided, the working directory
        and archive_dir are extended with the instance's workdir. This allows using descriptors to
        define data members that "inherit" the storage from their class. Might be removed.
    """

    def __init__(
        self,
        workdir: Path = None,
        archive_dir: Path = None,
        instance: StorableProtocol = None,
    ):
        self.workdir = Path(workdir if workdir else "").expanduser().absolute()
        self.archive_dir = (
            Path(archive_dir).expanduser().absolute() if archive_dir else None
        )
        self._member = None
        self.default_regex = []

        if instance is not None:
            self.workdir = self.workdir / instance.workdir
            if self.archive_dir is not None:
                self.archive_dir = self.archive_dir / instance.workdir

        if not self.workdir.exists():
            self.workdir.mkdir(parents=True)

    def glob(self, pattern: str) -> list[Path]:
        files = {}
        if self.archive_dir is not None:
            files.update(
                {
                    pth.relative_to(self.archive_dir): pth
                    for pth in self.archive_dir.glob(pattern)
                }
            )
        files.update(
            {pth.relative_to(self.workdir): pth for pth in self.workdir.glob(pattern)}
        )
        return list(files.values())

    def path(self, name: str) -> Path:
        if Path(name).is_absolute():
            raise ValueError(
                f"argument to Observation.file_path must be a relative path ({name})"
            )
        wp = self.workdir / name
        if self.archive_dir is not None:
            ap = self.archive_dir / name
            if not wp.exists() and ap.exists():
                return ap
        return wp

    def work_path(self, name: str) -> Path:
        name = Path(name)
        if self.archive_dir is not None and self.archive_dir in name.parents:
            name = name.relative_to(self.archive_dir)
        if name.is_absolute() and not (
            self.archive_dir in name.parents or self.workdir in name.parents
        ):
            raise ValueError(
                f"argument to Observation.file_path must either be a relative path ({name}) "
                f"or be in the archive or directories ({self.archive_dir}, {self.workdir})"
            )
        return self.workdir / name

    def archive(self, *regex: str):
        """
        Move files to archive location.

        Parameters
        ----------
        regex: list of str
            Optional. If not given, self.default_regex is used.
            Files matching any of the strings are arcived in a long-term location.
        """
        if self.archive_dir is None:
            raise RuntimeError("archive destination was not specified")

        regex = regex if regex else self.default_regex
        destination = self.archive_dir
        # logger.debug(f"{self} Archiving to {destination}:")
        for pattern in regex:
            for src in self.workdir.rglob(f"**/{pattern}"):
                # logger.debug(f"{self}   - {src}")
                dest = destination / src.relative_to(self.workdir)
                if not dest.parent.exists():
                    dest.parent.mkdir(exist_ok=True, parents=True)
                if src.is_dir():
                    if not dest.exists():
                        dest.mkdir(exist_ok=True, parents=True)
                    for src_2 in src.glob("**/*"):
                        if src_2.is_dir():
                            # Only files are copied. Parent directories are created automatically.
                            # empty directories are not copied
                            continue
                        dest_2 = destination / src_2.relative_to(self.workdir)
                        if not dest_2.parent.exists():
                            dest_2.parent.mkdir(exist_ok=True, parents=True)
                        try:
                            shutil.copy(src_2, dest_2)
                        except shutil.SameFileError:
                            # links to the same file are not copied
                            pass
                else:
                    try:
                        shutil.copy(src, dest)
                    except shutil.SameFileError:
                        pass

    def specialize(self, instance: StorableProtocol):
        """
        Return a Storage instance for a specific object.
        """
        return Storage(
            workdir=self.workdir,
            archive_dir=self.archive_dir,
            instance=instance,
        )

    def __repr__(self):
        workdir = self.workdir
        archive_dir = self.archive_dir
        return f"Storage({workdir=}, {archive_dir=})"


class StorableClass:
    """
    Implements the Storable protocol and manages a temporary work directory.

    It can be used in a derived class as follows:

        class PeriscopeDriftData(StorableClass):
            def __init__(self, obs):
            super().__init__(
                archive_dir=ARCHIVE_DIR,
                subdir=Path(f"{obs.obsid}")
            )

    The idea is that one has a top-level working directory and archive directory, and each
    instance gets its own subdirectory inside those. To some extent this is redundant with the
    instance argument of Storage.

    Parameters
    ----------
    workdir : str or Path, optional
        Path to the working directory. If None, a temporary directory is created.
    archive_dir : str or Path, optional
        Path to the archive directory. If None, no archiving is done.
    subdir : str, optional
        Subdirectory inside the workdir and archive_dir to use. Default is "".
    """

    def __init__(
        self,
        workdir=None,
        archive_dir=None,
        subdir="",
    ):
        self._clear = workdir is None
        self.tmp = tempfile.TemporaryDirectory() if workdir is None else None
        self.workdir = (
            Path(self.tmp.name if workdir is None else workdir).expanduser() / subdir
        )
        # checking if the workdir exists before creating is faster than passing the
        # exist_ok=True argument to mkdir, and there are thousands of observations
        # so this is a significant speedup
        if not self.workdir.exists():
            self.workdir.mkdir(parents=True)
        self.archive_dir = (
            (Path(archive_dir).expanduser() / subdir) if archive_dir else None
        )
        self.storage = Storage(
            workdir=self.workdir,
            archive_dir=self.archive_dir,
        )

    def __del__(self):
        if hasattr(self, "_clear") and self._clear:
            shutil.rmtree(self.workdir, ignore_errors=True)

    def file_glob(self, pattern):
        return self.storage.glob(pattern)

    def file_path(self, path):
        return self.storage.path(path)

    def archive(self, *regex: str):
        """
        Move files to archive location.

        Parameters
        ----------
        regex: list of str
            Optional. If not given, self.default_regex is used.
            Files matching any of the strings are arcived in a long-term location.
        """
        self.storage.archive(*regex)


class StoredResult:
    def __init__(
        self,
        *,
        func: Callable,
        name: str = None,
        fmt: str = None,
        storage: StorageProtocol = None,
        instance: StorableProtocol = None,
        subdir: Path | str = "",
        version: int = 0,
    ):
        """Cache `func`'s results under file names derived from its arguments.

        Parameters
        ----------
        version: int
            The version of what `func` returns. It is stored with every saved result,
            and a stored result from another version is a cache miss: it is recomputed
            and overwritten under the same file name. Bump it whenever a change alters
            the result -- otherwise results made by the old code keep being reused.
            Results stored before versions existed read back as version 0. A file
            written at a version above 0 cannot be read by astromon code from before
            versions existed (JSON and pickle results are wrapped with their version).
        """
        self.func = func  # the function being decorated
        self.version = int(version)
        self.name = (
            fullname(func) if name is None else name
        )  # the name is used to set the filename
        self.instance = instance
        self.fmt = "pickle" if fmt is None else fmt
        self.storage = Storage() if storage is None else storage
        self._member = None
        self.subdir = Path(subdir)

        if self.fmt not in FORMATS:
            raise ValueError(
                f"Unknown format {fmt}. Available formats: {list(FORMATS.keys())}"
            )

    def get_filename(self, *args, **kwargs) -> Path:
        if self.instance is not None:
            hsh = argument_hash(self.func, self.instance, *args, **kwargs)  # method
        else:
            hsh = argument_hash(self.func, *args, **kwargs)  # function
        # Use "__" rather than "::" to separate name from hash: macOS (APFS/HFS+)
        # disallows colons in filenames, which broke every cache write on macOS.
        filename = f"{self.name}__{hsh}"
        filename = FORMATS[self.fmt].sanitize_filename(filename)
        return self.subdir / filename

    def load(self, *args, **kwargs):
        filename = self.storage.path(self.get_filename(*args, **kwargs))
        if self.is_valid_result(*args, **kwargs) and filename.exists():
            result, version = FORMATS[self.fmt].load_versioned(filename)
            if version == self.version:
                return result
        return None

    def run(self, *args, **kwargs) -> Any:
        if self.instance is not None:
            result = self.func(self.instance, *args, **kwargs)  # method
        else:
            result = self.func(*args, **kwargs)  # function
        return result

    def save(self, result: Any, *args, **kwargs):
        filename = self.storage.work_path(self.get_filename(*args, **kwargs))
        FORMATS[self.fmt].save(result, filename, force=True, version=self.version)
        invalid = self._invalidate_filename(filename)
        if invalid.exists():
            invalid.unlink()

    def is_valid_result(self, *args, **kwargs) -> bool:
        filename = self.get_filename(*args, **kwargs)
        invalid = self._invalidate_filename(filename)
        return not invalid.exists()

    def invalidate_result(self, *args, **kwargs):
        filename = self.get_filename(*args, **kwargs)
        invalid = self._invalidate_filename(filename)
        if not invalid.parent.exists():
            invalid.parent.mkdir(parents=True, exist_ok=True)
        with open(invalid, "w"):
            pass

    def invalidate_all_results(self):
        for filename in self.storage.glob(f"{self.subdir / self.name}__*"):
            invalid = self._invalidate_filename(filename)
            if not invalid.parent.exists():
                invalid.parent.mkdir(parents=True, exist_ok=True)
            with open(invalid, "w"):
                pass

    def files(self) -> list[Path]:
        return list(self.storage.glob(f"{self.subdir}/{self.name}__*"))

    def _invalidate_filename(self, filename):
        filename = self.storage.work_path(filename)
        return filename.parent / f".{filename.name}.invalidated"

    def __call__(self, *args, **kwargs) -> Any:
        # check if there is a previous result and return it it exists
        result = self.load(*args, **kwargs)
        if result is None:
            # if there is no result, run the function and save the result
            result = self.run(*args, **kwargs)
            self.save(result, *args, **kwargs)
        return result

    def specialize(self, obj):
        """
        Return a StoredResult instance for a specific object.
        """
        return StoredResult(
            func=self.func,
            name=self.name,
            fmt=self.fmt,
            storage=obj.storage,
            subdir=self.subdir,
            instance=obj,
            version=self.version,
        )


class StoredResultWrapper(StoredResult):
    """
    A wrapper for StoredResult that allows it to be used as a descriptor.
    """

    def __get__(self, obj, objtype_):
        # This method is called when the StoredResultWrapper is an attribute of an object.
        # It returns a new instance of StoredResult with the object as the "instance".
        return self.specialize(obj)


def stored_result(*args, **kwargs) -> Callable:
    """
    Decorator to store the result of a function or method call.

        Usage:
        @store_result(name="my_function", fmt="json")
        def my_function(...):
            ...

    Pass ``version=N`` and bump it whenever a change alters what the function returns,
    so results cached by the old code are recomputed instead of reused (see
    :class:`StoredResult`). The version covers the code only: a cached function that
    reads a task's output must also notice when that task reruns, as
    ``Observation.get_sources`` does with the task run it stamps into each table.
    """

    if len(args) == 1 and callable(args[0]) and len(kwargs) == 0:
        # If the first argument is a callable, we are decorating a function
        return StoredResultWrapper(func=args[0])

    storage = kwargs.pop("storage", None)
    # Support both @stored_result("name", ...) and @stored_result(name="name", ...)
    name = args[0] if args and isinstance(args[0], str) else kwargs.pop("name", None)
    fmt = kwargs.pop("fmt", None)
    subdir = kwargs.pop("subdir", "")
    version = kwargs.pop("version", 0)

    def decorator(func: Callable) -> StoredResult:
        return StoredResultWrapper(
            func=func,
            name=name,
            fmt=fmt,
            storage=storage,
            subdir=subdir,
            version=version,
        )

    return decorator


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


def is_writable(dirpath):
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


def fullname(o):
    module = o.__class__.__module__
    if module == "builtins":
        return o.__qualname__  # avoid outputs like 'builtins.str'
    return f"{module}.{o.__qualname__}"


#: Where a non-zero result version is stored: a FITS header keyword for tables, a
#: wrapper key for JSON and pickle. Version 0 stores nothing, so a version-0 file is
#: exactly what was written before versions existed.
TABLE_VERSION_KEY = "RESVER"
WRAPPED_VERSION_KEY = "__stored_result_version__"


def _wrap(result: Any, version: int) -> Any:
    return {WRAPPED_VERSION_KEY: version, "value": result} if version else result


def _unwrap(stored: Any) -> tuple[Any, int]:
    if isinstance(stored, dict) and WRAPPED_VERSION_KEY in stored:
        return stored["value"], int(stored[WRAPPED_VERSION_KEY])
    return stored, 0


class TableCache:
    @staticmethod
    def save(table, filename, force=False, version=0):
        if filename.exists() and not force:
            return
        if not filename.parent.exists():
            filename.parent.mkdir(parents=True, exist_ok=True)
        # a shallow copy, so the caller's table keeps its meta as it was
        table = (
            table.copy(copy_data=False) if isinstance(table, Table) else Table(table)
        )
        if "description" in table.meta:
            table.meta["descript"] = table.meta["description"]
            del table.meta["description"]
        if version:
            table.meta[TABLE_VERSION_KEY] = int(version)
        table.write(filename, overwrite=force)

    @staticmethod
    def load_versioned(filename):
        table = Table.read(filename)
        table.convert_bytestring_to_unicode()
        version = int(table.meta.pop(TABLE_VERSION_KEY, 0))
        return table, version

    @staticmethod
    def load(filename):
        if filename.exists():
            return TableCache.load_versioned(filename)[0]

    @staticmethod
    def sanitize_filename(filename):
        filename = str(filename)
        if not filename.endswith((".fits", ".fits.gz")):
            filename = Path(str(filename) + ".fits.gz")
        return Path(filename)


class PickleCache:
    @staticmethod
    def save(result, filename, force=False, version=0):
        if filename.exists() and not force:
            return
        if not filename.parent.exists():
            filename.parent.mkdir(parents=True, exist_ok=True)
        with open(filename, "wb") as fh:
            pickle.dump(_wrap(result, version), fh)

    @staticmethod
    def load_versioned(filename):
        with open(filename, "rb") as fh:
            return _unwrap(pickle.load(fh))

    @staticmethod
    def load(filename):
        if filename.exists():
            return PickleCache.load_versioned(filename)[0]

    @staticmethod
    def sanitize_filename(filename):
        filename = str(filename)
        if not (filename.endswith((".pkl", ".pickle"))):
            filename = filename + ".pkl"
        return Path(filename)


class JsonCache:
    @staticmethod
    def save(result, filename, force=False, version=0):
        if filename.exists() and not force:
            return
        if not filename.parent.exists():
            filename.parent.mkdir(parents=True, exist_ok=True)
        with open(filename, "w") as fh:
            json.dump(_wrap(result, version), fh)

    @staticmethod
    def load_versioned(filename):
        with open(filename, "r") as fh:
            return _unwrap(json.load(fh))

    @staticmethod
    def load(filename):
        if filename.exists():
            return JsonCache.load_versioned(filename)[0]

    @staticmethod
    def sanitize_filename(filename):
        filename = str(filename)
        if not (filename.endswith(".json")):
            filename = filename + ".json"
        return Path(filename)


FORMATS = {
    "table": TableCache,
    "json": JsonCache,
    "pickle": PickleCache,
}
