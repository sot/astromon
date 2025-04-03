import os
import tempfile
from pprint import pprint
from unittest.mock import Mock

import pytest

from astromon import observation, stored_result, task


class MockTaskManager(task.TaskManager):
    """
    A mock task manager that has a Mock inner_function that can be called by tasks."""

    def __init__(self):
        super().__init__()
        self.inner_function = Mock()


def _call_args_list(mock):
    return [c.args for c in mock.call_args_list]


def get_obs(obsid):
    workdir = tempfile.TemporaryDirectory()
    archive_dir = tempfile.TemporaryDirectory()

    return observation.Observation(
        obsid, workdir=workdir.name, archive_dir=archive_dir.name
    )


@pytest.fixture()
def test_obs():
    yield get_obs("8007")


@pytest.fixture()
def test_pipeline():
    manager = MockTaskManager()

    @manager.task(
        name="one",
        outputs={
            "one": "test/one_1.txt",
            "two": "test/one_2.txt",  # this one is not created
            "three": "test/one_3.txt",
        },
    )
    def one(obs, inputs, outputs):
        manager.inner_function("one", obs.obsid)
        with open(outputs["one"], "w") as f:
            f.write("This is test file one for task one.\n")
        with open(outputs["three"], "w") as f:
            f.write("This is test file three for task one.\n")

    @manager.task(
        name="two",
        inputs={
            "one": "test/one_1.txt",
        },
        outputs={
            "two": "test/two_{obsid}_1.txt",
        },
    )
    def two(obs, inputs, outputs):
        manager.inner_function("two", obs.obsid)
        with open(inputs["one"], "r") as f:
            content = f.read()
        with open(outputs["two"], "w") as f:
            f.write(
                f"This is test file one for task two, reading from {inputs['one']}:\n{content}"
            )

    @manager.task(
        name="three",
        inputs={
            "one": "test/one_1.txt",
            # "two": "test/one_2.txt",
        },
        outputs={
            "one": "test/three_1.txt",
            "two": "test/three_2.txt",
        },
    )
    def three(obs, inputs, outputs):
        manager.inner_function("three", obs.obsid)
        with open(outputs["two"], "w") as f:
            f.write("This is test file two for task three.\n")

    @manager.task(
        name="four",
        inputs={
            "one": "test/two_{obsid}_1.txt",
            "two": "test/three_2.txt",
        },
        outputs={
            "one": "test/four_{obsid}_1.txt",
            "two": "test/four_{obsid}_2.txt",  # this one is not created
        },
    )
    def four(obs, inputs, outputs):
        manager.inner_function("four", obs.obsid)
        with open(inputs["one"], "r") as f:
            content = f.read()
        with open(outputs["one"], "w") as f:
            f.write(
                f"This is test file one for task two, reading from {inputs['one']}:\n{content}"
            )

    @manager.task(
        name="five",
        outputs={
            "two": "test/five_{band}_1.txt",
        },
        variables={
            "band": lambda obs: "wide" if obs.is_hrc else "broad",
        },
    )
    def five(obs, inputs, outputs):
        manager.inner_function("five", obs.obsid)
        with open(outputs["two"], "w") as f:
            f.write("This is test file one for task five")

    @manager.task(
        name="eight",
        outputs={
            "one": "test/eight.txt",  # this one is not created
        },
    )
    def eight(obs, inputs, outputs):
        manager.inner_function("eight", obs.obsid)

    @manager.task(
        name="six",
        inputs={
            "one": "test/four_{obsid}_1.txt",
        },
        optional_inputs={
            "two": "test/three_2.txt",
            "three": "test/eight.txt",  # this one is not created
        },
        outputs={
            "one": "test/six_{obsid}_1.txt",
        },
    )
    def six(obs, inputs, outputs):
        manager.inner_function("six", obs.obsid)
        with open(inputs["one"], "r") as f:
            content = f.read()
        with open(outputs["one"], "w") as f:
            f.write(
                f"This is test file one for task six, reading from {inputs['one']}:\n{content}"
            )

    @manager.task(
        name="seven",
        inputs={
            "one": "test/four_{obsid}_1.txt",
        },
        outputs={
            "one": "test/seven_{obsid}_1.txt",
        },
    )
    def seven(obs, inputs, outputs):
        manager.inner_function("seven", obs.obsid)
        with open(inputs["one"], "r") as f:
            content = f.read()
        with open(outputs["one"], "w") as f:
            f.write(
                f"This is test file one for task seven, reading from {inputs['one']}:\n{content}"
            )

    return manager


def test_a_precondition_dependency_list(test_pipeline):
    pprint(list(test_pipeline.task_graph.edges()))
    assert set(test_pipeline.task_graph.edges()) == {
        ("one", "two"),
        ("one", "three"),
        ("two", "four"),
        ("three", "four"),
        ("three", "six"),
        ("four", "six"),
        ("four", "seven"),
        ("eight", "six"),
    }, "Precondition for tests failed"


def test_dependency_check(test_pipeline):
    # test that the input files are checked
    obs = get_obs("8007")
    one = test_pipeline.tasks["one"]
    two = test_pipeline.tasks["two"]
    with pytest.raises(FileNotFoundError):
        two(obs)
    one(obs)
    two(obs)


def test_task_should_run_one(test_pipeline):
    # this function checks when a function should run
    # - if the cache is invalidated
    # - if the output files are missing (except ones that did not exist when the task was run)
    obs = get_obs("8007")
    obs2 = get_obs("8008")
    one = test_pipeline.tasks["one"]

    assert test_pipeline.tasks["one"].should_run(obs), "Task one should run initially"
    one(obs)
    assert not test_pipeline.tasks["one"].should_run(obs), (
        "Task one should not run after execution."
    )
    assert os.path.exists(obs.file_path("test/one_1.txt")), (
        "Output file one_1.txt should exist after task one execution."
    )
    assert not os.path.exists(obs.file_path("test/one_2.txt")), (
        "Output file one_2.txt should not exist after task one execution."
    )
    assert not test_pipeline.tasks["one"].should_run(obs, {"test/one_2.txt"}), (
        "Task one should not run after execution if it is a known missing file (1)."
    )
    assert not test_pipeline.tasks["one"].should_run(obs, {"test/one_3.txt"}), (
        "Task one should not run after execution."
    )

    # we now remove the output file to simulate a missing file
    obs.file_path("test/one_3.txt").unlink(missing_ok=True)
    assert not test_pipeline.tasks["one"].should_run(obs, {"test/one_2.txt"}), (
        "Task one should not run after execution if it is a known missing file (2)."
    )
    assert test_pipeline.tasks["one"].should_run(obs), (
        "Task one should run after execution if an output file is missing, even if cache is valid."
    )

    one.invalidate_result(obs)
    assert test_pipeline.tasks["one"].should_run(obs), (
        "Task one should run after cache clear."
    )
    one(obs)
    assert not test_pipeline.tasks["one"].should_run(obs), (
        "Task one should not run after cache clear."
    )
    one.invalidate_result(obs2)
    assert not test_pipeline.tasks["one"].should_run(obs), (
        "Task one should not run after cache clear on a different observation."
    )


def test_should_run_four(test_pipeline):
    obs = get_obs("8007")
    four = test_pipeline.tasks["four"]

    assert test_pipeline.inner_function.call_count == 0
    assert four.should_run(obs), "Task four should run initially"

    # we run task four dependencies by hand so we don't get errors in dependency resolution
    test_pipeline.tasks["one"](obs)
    test_pipeline.tasks["two"](obs)
    test_pipeline.tasks["three"](obs)

    assert four.should_run(obs), "Task four should run initially"
    four(obs)
    assert not four.should_run(obs), "Task four should not run again"
    four.invalidate_result(obs)
    assert four.should_run(obs), "Task four should run again after cache clear"
    four(obs)
    obs.file_path("test/three_2.txt").unlink(missing_ok=True)
    assert not four.should_run(obs), (
        "Task four should run again if three_2.txt is missing (missing inputs do not invalidate the result)"
    )
    four(obs)
    assert not four.should_run(obs), "Task four should not run again (2)"


def test_sequence_one(test_pipeline, test_obs):
    # test that a task should run once and only once
    test_pipeline.run_task(test_obs, "one")
    assert _call_args_list(test_pipeline.inner_function) == [("one", "8007")], (
        "Task one should have run"
    )
    test_pipeline.run_task(test_obs, "one")
    assert _call_args_list(test_pipeline.inner_function) == [("one", "8007")], (
        "Task one should run only once"
    )


def test_sequence_two(test_pipeline, test_obs):
    # task two depends on task one, so both should run
    test_pipeline.run_task(test_obs, "two")
    assert _call_args_list(test_pipeline.inner_function) == [
        ("one", "8007"),
        ("two", "8007"),
    ], "Task one and two should be run"


def test_sequence_five(test_pipeline, test_obs):
    # task five does not depend on any other task (this test a disconnected dependency graph)
    test_pipeline.run_task(test_obs, "five")
    assert _call_args_list(test_pipeline.inner_function) == [("five", "8007")], (
        "Task five should run"
    )


def test_sequence_four(test_pipeline, test_obs):
    # checking a diamon dependency chain
    test_obs = get_obs("8007")
    test_pipeline.run_task(test_obs, "four")
    # the task dependency graph is a diamon, tasks two and three are in the same level
    # so they can be run in any order
    possible_stacks = [
        [("one", "8007"), ("two", "8007"), ("three", "8007"), ("four", "8007")],
        [("one", "8007"), ("three", "8007"), ("two", "8007"), ("four", "8007")],
    ]
    assert _call_args_list(test_pipeline.inner_function) in possible_stacks, (
        "Task one, two, three and four should be run"
    )


def test_invalidate_result_output_1(test_pipeline, test_obs):
    # This tests checks that removing an output file invalidates the result.
    test_pipeline.run_task(test_obs, "four")
    test_pipeline.inner_function.reset_mock()

    # removing an output file invalidates the result, and the task (and nothing else) should be run
    test_obs.file_path(f"test/four_{test_obs.obsid}_1.txt").unlink(missing_ok=True)
    test_pipeline.run_task(test_obs, "four")
    assert _call_args_list(test_pipeline.inner_function) == [("four", "8007")], (
        "Task four should be run again due to missing file four_{obsid}_1.txt"
    )


def test_invalidate_result_output_2(test_pipeline, test_obs):
    # this test checks that a missing output file that was not created when the task ran
    # does not invalidate the result
    test_pipeline.run_task(test_obs, "four")
    test_pipeline.inner_function.reset_mock()

    # removing an output file that was not created does not invalidate the result
    test_obs.file_path(f"test/four_{test_obs.obsid}_2.txt").unlink(missing_ok=True)
    test_pipeline.run_task(test_obs, "four")
    assert _call_args_list(test_pipeline.inner_function) == [], (
        "Task four should not run again due to missing file four_{obsid}_2.txt"
    )


def test_invalidate_result(test_pipeline, test_obs):
    # This test checks the behaviour when a task is explicitly invalidated.
    test_pipeline.run_task(test_obs, "four")
    test_pipeline.inner_function.reset_mock()

    # explicitly invalidating the result, the task (and nothing else) should be run
    test_pipeline.tasks["four"].invalidate_result(test_obs)
    test_pipeline.run_task(test_obs, "four")
    assert _call_args_list(test_pipeline.inner_function) == [("four", "8007")], (
        "Task four should be run again due to invalid cache"
    )


def test_invalidate_result_input(test_pipeline, test_obs):
    # This test checks that removing the INPUT of a task does not invalidate the result.
    test_pipeline.run_task(test_obs, "four")
    test_pipeline.inner_function.reset_mock()

    # removing an input file does not invalidate the result
    test_obs.file_path("test/three_2.txt").unlink(missing_ok=True)
    test_pipeline.run_task(test_obs, "four")
    assert test_pipeline.inner_function.call_count == 0, (
        "Task four should not be run again because missing inputs do not invalidate the result"
    )


def test_invalidate_result_dependency_input(test_pipeline, test_obs):
    # This test invalidates a task (four), which STRICTLY depends on the output of another (three).
    # That means that task three should run if the file is missing.
    test_pipeline.run_task(test_obs, "four")
    test_pipeline.inner_function.reset_mock()

    # removing an output file invalidates the result,
    test_obs.file_path(f"test/four_{test_obs.obsid}_1.txt").unlink(missing_ok=True)
    # and removing an input from a dependency causes the dependency to be run as well
    test_obs.file_path("test/three_2.txt").unlink(
        missing_ok=True
    )  # required by task four
    test_pipeline.run_task(test_obs, "four")
    assert _call_args_list(test_pipeline.inner_function) == [
        ("three", "8007"),
        ("four", "8007"),
    ], (
        f"Task three and four should be run again due to missing files three_2.txt and four_{test_obs.obsid}_1.txt"
    )


def test_invalidate_result_dependency_input_2(test_pipeline, test_obs):
    # This test invalidates a task (six), which OPTIONALLY depends on the output of another (three).
    # That means that task three should run if the file is missing.
    test_pipeline.run_task(test_obs, "six")
    test_pipeline.inner_function.reset_mock()

    # explicitly invalidating the result of task six
    test_pipeline.tasks["six"].invalidate_result(test_obs)
    # and removing an optional input from a dependency causes the dependency to be run as well
    test_obs.file_path("test/three_2.txt").unlink(
        missing_ok=True
    )  # optional input for task six
    test_pipeline.run_task(test_obs, "six")
    assert _call_args_list(test_pipeline.inner_function) == [
        ("three", "8007"),
        ("six", "8007"),
    ], (
        f"Task three and four should be run again due to missing files three_2.txt and four_{test_obs.obsid}_1.txt"
    )


def test_invalidate_result_dependency_input_3(test_pipeline, test_obs):
    # This test invalidates task six, which OPTIONALLY depends on one output of task eight.
    # That means that task eight should run if the file is missing, unless we know that task eight
    # has run before and the file was not produced.
    test_pipeline.run_task(test_obs, "six")
    assert ("eight", "8007") in _call_args_list(test_pipeline.inner_function), (
        "Task eight should have run initially because it produces an optional input for task six"
    )
    test_pipeline.inner_function.reset_mock()

    # explicitly invalidating the result of task six
    test_pipeline.tasks["six"].invalidate_result(test_obs)
    # and removing a non-existent optional input from a dependency
    # does not cause the dependency to be run
    test_obs.file_path("test/eight.txt").unlink(missing_ok=True)
    test_pipeline.run_task(test_obs, "six")
    assert _call_args_list(test_pipeline.inner_function) == [("six", "8007")], (
        f"Task three and four should be run again due to missing files three_2.txt and four_{test_obs.obsid}_1.txt"
    )


def test_invalidate_dependency(test_pipeline, test_obs):
    # this test invalidates one task (three), which has the side-effect of invalidating a dependent
    # task (four). It then runs both tasks separately

    test_pipeline.run_task(test_obs, "four")
    test_pipeline.inner_function.reset_mock()

    # invalidating the result of task three should cause four to be run again,
    # even though task four's cache is still "valid"
    test_pipeline.tasks["three"].invalidate_result(test_obs)
    test_pipeline.run_task(test_obs, "four")
    assert _call_args_list(test_pipeline.inner_function) == [
        ("three", "8007"),
        ("four", "8007"),
    ], (
        "Task three and four should be run again because the cache of task three was cleared"
    )


def test_invalidate_dependency_file(test_pipeline, test_obs):
    # this test invalidates one task (four) by removing its output,
    # it also removes an output of task three, but task four does not depend on it
    # so task three should not be run again

    test_pipeline.run_task(test_obs, "four")
    test_pipeline.inner_function.reset_mock()

    # Task four depends on three, but only through the output file three_2.txt,
    # removing the other output of task three should not cause task three to be run when re-running four
    test_obs.file_path("test/three_1.txt").unlink(missing_ok=True)
    test_obs.file_path(f"test/four_{test_obs.obsid}_1.txt").unlink(missing_ok=True)
    test_pipeline.run_task(test_obs, "four")
    assert _call_args_list(test_pipeline.inner_function) == [("four", "8007")], (
        f"Task four should be run again because missing four_{test_obs.obsid}_1.txt and task three should not run even though three_1.txt is missing"
    )


def test_invalidate_dependency_two_steps(test_pipeline, test_obs):
    # this test invalidates one task (three), which has the side-effect of invalidating a dependent
    # task (four). It then runs both tasks separately
    test_pipeline.run_task(test_obs, "four")
    test_pipeline.inner_function.reset_mock()

    # invalidate cache in the middle, only run the intermediate task, then run the last task
    test_pipeline.tasks["three"].invalidate_result(test_obs)
    test_pipeline.run_task(test_obs, "three")  # three is re-run here, four is not
    assert _call_args_list(test_pipeline.inner_function) == [("three", "8007")], (
        "Task three should run again because its cache was cleared"
    )

    # four should run here even though its result has not been explicitly invalidated
    test_pipeline.run_task(test_obs, "four")
    assert _call_args_list(test_pipeline.inner_function) == [
        ("three", "8007"),
        ("four", "8007"),
    ], "Task four should run again because the cache of task three was cleared"


def test_invalidate_archive(test_pipeline, test_obs):
    # running "four" to begin with
    test_pipeline.run_task(test_obs, "four")
    test_pipeline.inner_function.reset_mock()

    # archiving all but task three
    test_obs.archive("test/one*", "test/two*", "test/four*", "cache/*")

    # running "four" again should not run the task, because the result is still valid
    test_pipeline.run_task(test_obs, "four")
    assert _call_args_list(test_pipeline.inner_function) == [], "Task already ran"

    # create an observation with a different work directory but the same archive directory
    workdir = tempfile.TemporaryDirectory()
    obs_2 = observation.Observation(
        test_obs.obsid,
        workdir=workdir.name,
        archive_dir=test_obs.archive_dir.parent.parent,
    )

    os.chmod(obs_2.archive_dir, 0o500)  # ensure the archive is NOT writable
    with pytest.raises(PermissionError):
        open(obs_2.archive_dir / "fail.txt", "w").close()

    # task "four" should not run on the new observation,
    # because the result is still valid in the archive directory
    test_pipeline.run_task(obs_2, "four")
    stack = _call_args_list(test_pipeline.inner_function)
    assert stack == [], "Task already ran (archived)"

    test_pipeline.tasks["four"].invalidate_result(obs_2)

    # after invalidating the result, the task should run again
    # and because we did not archive task three, it should run as well
    test_pipeline.run_task(obs_2, "four")
    stack = _call_args_list(test_pipeline.inner_function)
    assert stack == [("three", "8007"), ("four", "8007")], (
        "Tasks three and four should run again"
    )


def test_task_return_value(test_pipeline):
    # check the possible return values: None, and tuples of size 1 and 2

    @test_pipeline.task()
    def none(obs, inputs, outputs):
        pass

    @test_pipeline.task()
    def one(obs, inputs, outputs):
        return task.ReturnCode.OK

    @test_pipeline.task()
    def two(obs, inputs, outputs):
        return task.ReturnCode.SKIP, "skipped"

    obs = get_obs("8007")

    rv = none(obs)
    assert isinstance(rv, task.ReturnValue), (
        "Return value should be a ReturnValue object"
    )
    assert rv

    rv = two(obs)
    assert isinstance(rv, task.ReturnValue), (
        "Return value should be a ReturnValue object"
    )
    assert not rv
    assert rv.return_code == task.ReturnCode.SKIP
    assert rv.msg == "skipped"

    rv = one(obs)
    assert isinstance(rv, task.ReturnValue), (
        "Return value should be a ReturnValue object"
    )
    assert rv
    assert rv.return_code == task.ReturnCode.OK
    assert rv.msg == ""


def test_download(test_pipeline):
    # check that some files are downloaded from the archive
    @test_pipeline.task(
        name="dnld",
        download=["evt2", "fov"],
    )
    def dnls(obs, inputs, outputs):
        pass

    obs = get_obs("8007")
    dnls(obs)
    assert len(list(obs.workdir.glob("*/*evt2*"))) == 1, "missing event files"
    assert len(list(obs.workdir.glob("*/*fov*"))) == 1, "missing FOV files"


def test_repeated_task_name(test_pipeline):
    # check that defining a task with the same name raises an error
    with pytest.raises(ValueError):

        @test_pipeline.task(
            name="five",
        )
        def five(obs, inputs, outputs):
            pass


def test_dependencies():
    # Test the "dependencies" decorator with a simple task

    TASKS = task.TaskManager()
    STACK = []
    TMPDIR = tempfile.TemporaryDirectory()

    @TASKS.task(
        name="do",
        outputs={"do": "do_nothing.json"},
    )
    def do(obs, inputs=None, outputs=None):
        STACK.append(("do", obs.obsid))
        with open(outputs["do"], "w") as f:
            f.write("This is a placeholder task that does nothing.")

    class Data:
        def __init__(self):
            self.storage = stored_result.Storage(workdir=TMPDIR.name)
            self.obsid = "1"

        @property
        def workdir(self):
            """
            Returns the work directory for the current observation.
            """
            return self.storage.workdir

        def file_path(self, *args, **kwargs):
            """
            Returns the file path for the given arguments.
            """
            return self.storage.path(*args, **kwargs)

        def file_glob(self, *args, **kwargs):
            """
            Returns the file path for the given arguments.
            """
            return self.storage.glob(*args, **kwargs)

        def download(self, *args, **kwargs):
            """
            Downloads the file for the given arguments.
            """

        @TASKS.dependencies(required_files={"do": "do_{name}.json"})
        def method(self, name="nothing", apply_filter=True):
            return {"name": name, "apply_filter": apply_filter}

    data = Data()
    # check that the method runs and returns the expected result
    assert data.method() == {"name": "nothing", "apply_filter": True}
    # check that the result is the same with the default argument
    assert data.method(name="nothing") == {"name": "nothing", "apply_filter": True}
    # check that the "do" task ran once and only once
    assert STACK == [("do", "1")], "Task do should be run exactly once"

    # the method does not accept positional arguments
    with pytest.raises(RuntimeError, match="positional arguments"):
        data.method("nothing")

    # test that requiring a file that is not produced by any tasks raises an exception
    # (no task produces do_bla.json)
    with pytest.raises(FileNotFoundError, match="do_bla.json"):
        data.method(name="bla")  # fails, file not found
    assert STACK == [("do", "1")], "Task do should be run exactly once"
