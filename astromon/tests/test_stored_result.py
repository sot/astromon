import contextlib
import tempfile
import time
from pathlib import Path

import numpy as np
import pytest

from astromon.stored_result import Storage, StoredResult, stored_result


def get_sum(a, b):
    return a + b


def get_diff(a, b):
    return a - b


def get_log(a, base=10):
    return np.log(a) / np.log(base)


def get_random(seed=None):
    return np.random.rand()


class GetLog:
    def get_log_member(self, a, base=10):
        return get_log(a, base)


get_log_instance = GetLog()


@pytest.fixture
def temp_storage():
    with tempfile.TemporaryDirectory() as td:
        storage = Storage(workdir=td)
        yield storage


@pytest.fixture
def temp_archived_storage():
    with tempfile.TemporaryDirectory() as work, tempfile.TemporaryDirectory() as arch:
        storage = Storage(workdir=work, archive_dir=arch)
        yield storage


# this fixture is used automatically for all tests to avoid polluting the current working directory
@pytest.fixture(autouse=True)
def tmp_cwd():
    with tempfile.TemporaryDirectory() as td, contextlib.chdir(td):
        yield Path(td)


@pytest.fixture
def tmp_archive():
    with tempfile.TemporaryDirectory() as td, contextlib.chdir(td):
        yield Path(td)


def test_stored_result_constructor():
    # should raise TypeError
    with pytest.raises(TypeError):
        StoredResult()

    # should raise TypeError (no positional arguments)
    with pytest.raises(TypeError):
        StoredResult(get_random)

    sr1 = StoredResult(func=get_random)
    assert sr1.name == "get_random"

    sr2 = StoredResult(func=get_random, name="func.get_random")
    assert sr2.name == "func.get_random"

    sr1 = StoredResult(func=get_random, subdir="cache")
    assert sr1.name == "get_random"
    assert f"{sr1.subdir}" == "cache"

    # NOTE: this is not exactly how it's called
    sr3 = StoredResult(func=GetLog.get_log_member, instance=get_log_instance)
    assert sr3.name == "GetLog.get_log_member", (
        "Name should be derived from the instance method"
    )
    assert sr3.instance == get_log_instance, "Instance should be set correctly"

    # NOTE: this could have raised an exception as the "instance" does not have the method
    # with pytest.raises(AttributeError):
    #     StoredResult(func=get_log_instance.get_log, instance=get_log)


def test_format_filename(temp_storage):
    sr = StoredResult(func=get_sum, storage=temp_storage, fmt="json")
    assert sr.fmt == "json", "Format should be set to 'json'"
    assert sr.get_filename(1, 2).suffix == ".json", "Filename should have .json suffix"

    sr = StoredResult(func=get_diff, storage=temp_storage, fmt="pickle")
    assert sr.fmt == "pickle", "Format should be set to 'pickle'"
    assert sr.get_filename(1, 2).suffix == ".pkl", "Filename should have .pkl suffix"

    sr = StoredResult(func=get_log, storage=temp_storage, fmt="table")
    assert sr.fmt == "table", "Format should be set to 'table'"
    assert sr.get_filename(1, 2).name.endswith(".fits.gz"), (
        "Filename should have .fits.gz suffix"
    )


def test_filename(temp_storage):
    sr1 = StoredResult(func=get_log, storage=temp_storage)
    path_1_1 = sr1.get_filename(2)
    path_1_2 = sr1.get_filename(2, base=10)
    path_1_3 = sr1.get_filename(2, 10)

    assert isinstance(path_1_1, Path), "Filename should be a Path object"
    assert not path_1_1.is_absolute(), "Filename should not be absolute"
    assert path_1_1.name.startswith("get_log"), (
        "Filename should start with the function name"
    )

    assert path_1_1 == path_1_2, (
        "Filename should be the same for different keyword arguments"
    )
    assert path_1_2 == path_1_3, (
        "Filename should be the same regardless of keyword or positional arguments"
    )

    path_2 = sr1.get_filename(2, base=2)
    assert path_1_1 != path_2, "Filename should be different for different arguments"


def test_files():
    sr = StoredResult(func=get_sum)
    assert len(sr.files()) == 0, (
        "Files should be empty before the function is ever called"
    )
    sr(1, 2)
    sr(3, 4)
    assert len(sr.files()) == 2, (
        "Files should contain two files after the function is called"
    )
    for f in sr.files():
        # NOTE: this is different from the get_filename method, which returns relative path!
        assert f.is_absolute(), "File paths should be absolute"
        assert f.name.startswith("get_sum"), (
            "File name should start with the function name"
        )


def test_call():
    # the "call" method tries to load the result and calls the function and saves the result if needed
    sr1 = StoredResult(func=get_log)
    assert sr1.storage.path(sr1.get_filename(2)).exists() is False, (
        "Filename should not exist before the function is called"
    )
    assert sr1(2) == get_log(2), (
        "StoredResult should return the same value returned by the function"
    )
    # the following is True only because of the tmp_cwd fixture, which is used automatically
    # and changes the current working directory.
    assert sr1.storage.path(sr1.get_filename(2)).exists(), (
        "Filename should exist after the function is called"
    )


def test_run():
    # NOTE: this works. Is it the required behavior?
    # the "run" method calls the function directly (no save/load)
    sr1 = StoredResult(func=get_log)
    sr2 = StoredResult(func=GetLog.get_log_member, instance=get_log_instance)
    assert sr1.run(2) == get_log(2)
    assert sr1.run(2) == sr2.run(2)


def test_save_n_load():
    sr1 = StoredResult(func=get_log)
    assert sr1.storage.path(sr1.get_filename(2)).exists() is False, (
        "Filename should not exist before the function is called"
    )

    # NOTE: this one is currently returning None instead of raising an error
    # with pytest.raises(FileNotFoundError):
    #     sr1.load(2, base=10)

    sr1.save(125, 2, base=10)
    assert sr1.load(2, base=10) == 125
    assert sr1(2) == 125, "Function should return the correct value"

    sr1.invalidate_result(2, base=10)
    # NOTE: the following should raise some sort of error, because the result is invalid,
    # but currently it returns None, so it is not possible to know the difference between
    # and invalid result and a function returning None
    # with pytest.raises(RuntimeError):
    #     sr1.load(2, base=10)


@pytest.mark.parametrize(["subdir"], [("",), ("cache",)])
def test_invalidate(subdir):
    sr = StoredResult(func=get_random, subdir=subdir)
    sr.invalidate_result()  # invalidate a result that does not exist should be no-op

    value = sr()
    assert sr() == value, "StoredResult should return the same value as the first call"
    sr.invalidate_result()
    # NOTE: this fails because invalidate does not remove the file
    # assert len(sr.files()) == 0, "Files should contain no files after invalidated"
    value_2 = sr()
    assert value_2 != value, (
        "StoredResult should return a different value after invalidation"
    )
    # calling again in case the function ran last time and gave the right result, but wasn't cached
    assert sr() == value_2, (
        "StoredResult should return the first value returned after invalidation"
    )


@pytest.mark.parametrize(["subdir"], [("",), ("cache",)])
def test_invalidate_all(subdir):
    sr = StoredResult(func=get_random, subdir="")
    assert len(sr.files()) == 0, "Files should contain no files at the start"
    sr.invalidate_all_results()  # invalidate all results that do not exist should be no-op
    # sanity check

    value_1 = sr()
    value_2 = sr(2)
    assert value_1 != value_2, "StoredResult should return different values"
    assert sr() == value_1, (
        "StoredResult should return the same value as the first call"
    )
    assert sr(2) == value_2, (
        "StoredResult should return the same value as the second call"
    )
    assert len(sr.files()) == 2, (
        "Files should contain two files after the function is called"
    )
    sr.invalidate_all_results()
    # NOTE: this fails because invalidate does not remove the file
    # assert len(sr.files()) == 0, "Files should contain no files after invalidated"
    value_1_2 = sr()
    value_2_2 = sr(2)
    assert value_1_2 != value_1, (
        "StoredResult should return a different value after invalidation"
    )
    assert value_2_2 != value_2, (
        "StoredResult should return a different value after invalidation"
    )
    # calling again in case the function ran last time and gave the right result, but wasn't cached
    assert sr() == value_1_2, (
        "StoredResult should return the first value returned after invalidation"
    )
    assert sr(2) == value_2_2, (
        "StoredResult should return the second value returned after invalidation"
    )


def test_stored_result_method():
    class Data:
        def __init__(self, i):
            self.workdir = f"whatever_{i}"
            self.storage = Storage("test_dir_1", instance=self)

        @stored_result
        def method_1(self):
            pass

        @stored_result
        def method_2(self):
            pass

        @stored_result
        def my_function(self, a, b):
            """A simple function that adds two numbers."""
            time.sleep(2)
            return a - b

    data = Data(1)
    assert data.method_1.name != data.method_2.name
    assert data.method_1.name != data.my_function.name

    data_2 = Data(2)

    assert data.storage.workdir != data_2.storage.workdir, (
        "Workdirs should be different for different instances"
    )

    # directories are empty to begin with
    assert not data.storage.path(data.my_function.get_filename(1, 2)).exists()

    filename = data.storage.path(data.my_function.get_filename(1, 2))
    assert not filename.exists()
    start = time.time()
    data.my_function(1, 2)  # This will call the function and save the result
    stop = time.time()
    assert stop - start > 2
    assert filename.exists()
    start = time.time()
    data.my_function(1, 2)  # This will call the function and save the result
    stop = time.time()
    assert stop - start < 0.1
    assert filename.exists()

    # just checking that the filename portion of the path depends only on args and kwargs
    storage = Storage("test_dir_1")

    @stored_result(storage=storage)
    def my_function(a, b):
        """A simple function that adds two numbers."""
        return a - b

    assert str(data.my_function.get_filename(1, 2)) == str(
        my_function.get_filename(1, 2)
    ).replace("my_function", "Data.my_function")


def test_stored_result_function():
    storage = Storage("test_dir_1", "test_dir_2")

    @stored_result(storage=storage, name="my_function")
    def my_function(x, y):
        """
        A simple function that adds two numbers.
        """
        time.sleep(2)
        return x + y

    filename = storage.path(my_function.get_filename(1, 2))
    assert not filename.exists()
    start = time.time()
    my_function(1, 2)  # This will call the function and save the result
    stop = time.time()
    assert stop - start > 2
    assert filename.exists()
    assert len(my_function.files()) == 1, "Function files should contain one file"
    assert str(my_function.get_filename(1, 2)) == my_function.files()[0].name, (
        "Function files should contain the result file"
    )

    start = time.time()
    my_function(1, 2)  # This will call the function and save the result
    stop = time.time()
    assert stop - start < 0.1

    # archive and remove the file from the workdir
    storage.archive("my_function::*")
    filename.unlink(missing_ok=True)  # remove the file from the workdir

    assert len(my_function.files()) == 1, "Function files should contain one file"
    assert storage.archive_dir in my_function.files()[0].parents, (
        "Function files should be in the archive directory now"
    )
    assert str(my_function.get_filename(1, 2)) == my_function.files()[0].name, (
        "Function files should contain the result file"
    )

    start = time.time()
    my_function(1, 2)
    stop = time.time()
    assert stop - start < 0.1
    # the result is stored in the archive
    assert storage.path(my_function.get_filename(1, 2)).exists()
    assert storage.archive_dir in storage.path(my_function.get_filename(1, 2)).parents

    # invalidate the result
    my_function.invalidate_result(1, 2)

    # the result is still there (invalidating the result should not change the archive)
    assert storage.path(my_function.get_filename(1, 2)).exists()

    # now it should take time again
    start = time.time()
    my_function(1, 2)
    stop = time.time()
    assert stop - start > 2, "invalidate_result call did not invalidate result"

    start = time.time()
    my_function(1, 2)
    stop = time.time()
    assert stop - start < 0.1, "value not stored after invalidate_result call"

    # invalidate all results
    my_function.invalidate_all_results()

    # now it should take time again
    start = time.time()
    my_function(1, 2)
    stop = time.time()
    assert stop - start > 2, "invalidate_all_results call did not invalidate result"

    start = time.time()
    my_function(1, 2)  # This will call the function and save the result
    stop = time.time()
    assert stop - start < 0.1, "value not srored after invalidate_all_results call"


def test_stored_result_method_args():
    storage = Storage("test_dir_1")

    class Data:
        def __init__(self, i=0):
            self.workdir = f"whatever_{i}"
            self.storage = storage.specialize(instance=self)

        @stored_result(name="name", fmt="json", subdir="cache")
        def name(self, a, b):
            """A simple function that adds two numbers."""
            time.sleep(2)
            return (a, b)

    data = Data()
    dat = data.name("a", 1)

    assert dat == ("a", 1)
    assert len(data.name.files()) > 0, "Cache directory should contain the result file"


def test_storage_path():
    with tempfile.TemporaryDirectory() as td, contextlib.chdir(td):
        storage = Storage("test_dir_1", "archive_dir_1")

        path_w = storage.path("test_file.txt")  # should be in the workdir
        assert path_w.is_absolute(), "Path should be absolute"
        assert storage.workdir in path_w.parents, "Path should be in the workdir"

        with open(path_w, "w") as fh:
            fh.write("Hello, world!")
        storage.archive("test_file.txt")  # archive the file

        assert storage.workdir in storage.path("test_file.txt").parents, (
            "Path should be in workdir (takes precedence over archive)"
        )
        path_w.unlink()

        path_a = storage.path("test_file.txt")
        assert path_a.is_absolute(), "Path should be absolute"
        assert storage.archive_dir in path_a.parents, "Path should be in the archive"

        work_path = storage.work_path("test_file.txt")
        assert work_path == path_w
        work_path = storage.work_path(path_a)
        assert work_path == path_w


def test_storage():
    with tempfile.TemporaryDirectory() as td, contextlib.chdir(td):
        storage = Storage("test_dir_1")
        assert storage.workdir == Path("test_dir_1").absolute()
        assert storage.workdir.exists()
        assert storage.archive_dir is None

        with pytest.raises(RuntimeError):
            storage.archive()  # raises exception because the archive is not given

        # creating a storage with workdir and archive_dir
        storage = Storage("test_dir_1", "test_dir_2")
        assert storage.workdir == Path("test_dir_1").absolute()
        assert storage.archive_dir == Path("test_dir_2").absolute()
        filename = storage.path("test_1.txt")
        with open(filename, "w") as fh:
            fh.write("Hello, world!")

        # archiving the default (nothing is archive)
        storage.archive()
        assert not storage.archive_dir.exists()

        # archiving the file now
        storage.archive("test_*")
        assert storage.archive_dir.exists()

        # Create a new storage with same archive, different workdir
        storage = Storage("test_dir_3", "test_dir_2")
        assert storage.archive_dir.exists()
        # now the file should be in the archive
        filename = storage.path("test_1.txt")
        assert filename.exists()
        assert filename == storage.archive_dir / "test_1.txt"
        assert storage.archive_dir in filename.parents
        with open(storage.archive_dir / "test_1.txt", "r") as fh:
            assert fh.read() == "Hello, world!"

        # create a second file
        filename = storage.workdir / "test_2.txt"
        with open(filename, "w") as fh:
            fh.write("Hello, again!")

        # a glob should return both files (one in the archive, one in the workdir)
        assert storage.glob("test*") == [
            Path("test_dir_2/test_1.txt").absolute(),
            Path("test_dir_3/test_2.txt").absolute(),
        ]

        # overwriting the first file in the workdir
        with open(storage.workdir / "test_1.txt", "w") as fh:
            fh.write("Hello, world!")

        # now the glob should return both files in the workdir
        assert storage.glob("test*") == [
            Path("test_dir_3/test_1.txt").absolute(),
            Path("test_dir_3/test_2.txt").absolute(),
        ]

        class Data:
            storage = Storage("test_dir_1")

        data = Data()
        # with pytest.raises(AttributeError):
        #     _ = data.storage.workdir  # raises AttributeError: 'Data' object has no attribute 'workdir'

    with tempfile.TemporaryDirectory() as td, contextlib.chdir(td):
        # now check that workdir and archive_dir are set correctly in the instance
        class Data:
            def __init__(self):
                self.workdir = "whatever"
                self.storage = Storage("test_dir_1", instance=self)

        data = Data()
        assert data.storage.workdir == Path("test_dir_1/whatever").absolute()
        assert data.storage.archive_dir is None

        class Data:
            def __init__(self):
                self.workdir = "whatever"
                self.storage = Storage("test_dir_1", "test_dir_2", instance=self)

        data = Data()

        assert data.storage.workdir == Path("test_dir_1/whatever").absolute()
        assert data.storage.archive_dir == Path("test_dir_2/whatever").absolute()
