import dataclasses
import functools
import inspect
import json
import logging
import weakref
from collections import defaultdict
from enum import Enum
from pathlib import Path

import networkx as nx
from cxotime import CxoTime

logger = logging.getLogger("astromon")


class ReturnCode(Enum):
    """
    Enum to represent the return codes of a function.
    """

    OK = 0
    SKIP = 1
    WARNING = 2
    ERROR = 3
    INVALID = 1000  # for invalidating the result


@dataclasses.dataclass
class ReturnValue:
    """
    Class to encapsulate the return value of a function.

    Parameters
    ----------
    return_code: ReturnCode
        The return code of the function call. One of OK, SKIP, WARNING or ERROR.
    msg: str
        Optional message to be logged or displayed.
    """

    return_code: ReturnCode = ReturnCode.OK
    msg: str = ""
    actual_outputs: list = dataclasses.field(default_factory=list)
    inputs: dict = dataclasses.field(default_factory=dict)
    outputs: dict = dataclasses.field(default_factory=dict)
    start: CxoTime = dataclasses.field(default_factory=CxoTime.now)
    stop: CxoTime = dataclasses.field(default_factory=CxoTime.now)

    def __bool__(self):
        return self.return_code == ReturnCode.OK

    def to_dict(self):
        """
        Convert the ReturnValue to a JSON serializable dictionary.
        """
        return {
            "return_code": self.return_code.name,
            "msg": self.msg,
            "actual_outputs": self.actual_outputs,
            "inputs": self.inputs,
            "outputs": self.outputs,
            "start": self.start.iso,
            "stop": self.stop.iso,
        }

    @staticmethod
    def from_dict(data):
        """
        Populate the ReturnValue from a JSON serializable dictionary.
        """
        rv = ReturnValue()
        rv.return_code = ReturnCode[data["return_code"]]
        rv.msg = data.get("msg", "")
        rv.actual_outputs = data.get("actual_outputs", [])
        rv.inputs = data.get("inputs", {})
        rv.outputs = data.get("outputs", {})
        rv.start = CxoTime(data["start"])
        rv.stop = CxoTime(data["stop"])
        return rv


class Dependent:
    """
    Decorator for functions that depend on tasks.
    """

    def __init__(
        self,
        func,
        manager,
        instance=None,
        name=None,
        tasks=None,
        required_files=None,
        optional_files=None,
        download=None,
        variables=None,
    ):
        """
        Decorator to run tasks.

        Parameters
        ----------
        func: callable
            The function to be decorated.
        name: str
            Name of the task.
        tasks: list
            List of task names that this function depends on. An exception is raised if any of
            these tasks give an error.
        required_files: dict
            Dictionary of required input files. An exception is raised if any of these files are
            missing.
        optional_files: dict
            Dictionary of optional input files.
        download: list
            List of file categories to download before running the task.
        variables: dict
            Dictionary of variables to be used in the task. If a value is a function, it will be
            called without arguments.
        """
        self._instance = instance
        self._manager = manager
        self.func = func
        self.name = _fullname(func) if name is None else name
        self._tasks = tasks if tasks is not None else []
        self._required_files = required_files if required_files is not None else {}
        self._optional_files = optional_files if optional_files is not None else {}
        self._download = download if download is not None else []
        self._variables = variables if variables is not None else {}
        functools.update_wrapper(self, func)

    def __call__(self, *args, **kwargs):
        """
        Call the function
        """
        if self._instance is None:
            obs = args[0]
            args = args[1:]
        else:
            obs = self._instance

        if len(args) != 0:
            raise RuntimeError(
                "Dependent functions should not take positional arguments."
            )

        if self._download:
            obs.download(self._download)

        params = self.get_parameters(obs, **kwargs)

        requested_files = list(
            set(params["required_files"].values())
            | set(params["optional_files"].values())
        )

        rv = self._manager.run_tasks(obs=obs, requested_files=requested_files)

        errors = {
            name: value
            for name, value in rv.items()
            if value.return_code.value >= ReturnCode.ERROR.value
        }
        if errors:
            msg = ", ".join(f"{name} {value.msg}" for name, value in errors.items())
            raise RuntimeError(f"{self.name} failed. Dependency tasks failed: {msg}")

        missing = {
            name: value
            for name, value in params["required_files"].items()
            if not obs.file_path(value).exists()
        }
        if missing:
            msg = ", ".join(f"{value}" for value in missing.values())
            raise FileNotFoundError(f"{self.name} failed. Missing files: {msg}")

        return self.func(obs, **kwargs)

    def get_parameters(self, obs, **kwargs):
        """
        Get the task parameters for the observation, including variable interpolation.
        """
        signature = inspect.signature(self.func)
        arguments = signature.bind(obs, **kwargs)
        arguments.apply_defaults()

        variables = {"obsid": obs.obsid}
        variables.update(
            {
                key: value(obs) if callable(value) else value
                for key, value in self._variables.items()
            }
        )
        variables.update(arguments.arguments)

        parameters = {
            "required_files": {
                key: value.format(**variables)
                for key, value in self._required_files.items()
            },
            "optional_files": {
                key: value.format(**variables)
                for key, value in self._optional_files.items()
            },
            "variables": variables,
            "download": self._download,
        }
        return parameters


class Task:
    """
    Decorator to run tasks.
    """

    def __init__(
        self,
        func,
        manager,
        name=None,
        inputs=None,
        optional_inputs=None,
        outputs=None,
        download=None,
        variables=None,
    ):
        """
        Decorator to run tasks.

        Parameters
        ----------
        func: callable
            The function to be decorated.
        name: str
            Name of the task.
        inputs: dict
            Dictionary of input files.
        outputs: dict
            Dictionary of output files.
        download: list
            List of file categories to download before running the task.
        variables: dict
            Dictionary of variables to be used in the task. If a value is a function, it will be
            called without arguments.
        """
        self._manager = weakref.ref(manager)
        self.subdir = Path("cache")
        self.func = func
        self.name = _fullname(func) if name is None else name
        self._inputs = inputs if inputs is not None else {}
        self._optional_inputs = optional_inputs if optional_inputs is not None else {}
        self._outputs = outputs if outputs is not None else {}
        self._download = download if download is not None else []
        self._variables = variables if variables is not None else {}
        functools.update_wrapper(self, func)

    def __call__(self, obs, requested=None):
        """
        Same as Task.run
        """
        return self.run(obs, requested=requested)

    def get_parameters(self, obs):
        """
        Get the task parameters for the observation, including variable interpolation.
        """
        variables = {"obsid": obs.obsid}
        variables.update(
            {
                key: value(obs) if callable(value) else value
                for key, value in self._variables.items()
            }
        )
        parameters = {
            "inputs": {
                key: value.format(**variables) for key, value in self._inputs.items()
            },
            "optional_inputs": {
                key: value.format(**variables)
                for key, value in self._optional_inputs.items()
            },
            "outputs": {
                key: value.format(**variables) for key, value in self._outputs.items()
            },
            "variables": variables,
            "download": self._download,
        }
        return parameters

    def get_filename(self, obs):
        return obs.file_path(self.subdir / f"task_{self.name}.json")

    def invalidate_result(self, obs, follow_dependents=True):
        """
        Invalidate the result of the task.

        Parameters
        ----------
        obs: Observation
            The observation object.
        follow_dependents: bool
            If True, invalidate the result of all dependent tasks as well.
            If False, only invalidate the result of this task.
        """
        manager = self._manager()
        if manager is None:
            # we need the manager because we need to know the dependencies
            raise RuntimeError(
                "Task manager is not available. Cannot invalidate result."
            )

        # the cache file can be in the archive or the workdir, so we need to check both
        # the cache value is not removed, but set to invalid. This is because we can't always
        # write to the archive directory, but we might still override the cache.
        rv = self.get_result(obs)
        # if rv is None, there is nothing to clear
        if rv is not None:
            rv.return_code = ReturnCode.INVALID
            self.set_result(obs, rv)

        logger.debug(f"{obs} Cleared result of task {self.name} in {obs.workdir}.")

        if follow_dependents:
            for task in manager.get_dependents(self.name).values():
                task.invalidate_result(obs, follow_dependents=False)

    def get_result(self, obs):
        """
        Get the stored result of the task.

        Parameters
        ----------
        obs: Observation
            The observation object.

        Returns
        -------
        ReturnValue or None
            The stored result of the task, or None if the cache is corrupted or does not exist.
        """
        # reading from either the archive or work directories
        cache_file = obs.file_path(self.subdir / f"task_{self.name}.json")
        if cache_file.exists():
            with open(cache_file, "r") as fh:
                try:
                    return ReturnValue.from_dict(json.load(fh))
                except json.JSONDecodeError as exc:
                    logger.warning(
                        f"Cache file {cache_file} is corrupted "
                        f"({exc}). Task {self.name} must run again."
                    )
                    return None

    def set_result(self, obs, return_value):
        """
        Set the stored result of the task.

        Parameters
        ----------
        obs: Observation
            The observation object.
        return_value: ReturnValue
            The return value to be stored.
        """
        # never writing to the archive directory
        cache_file = obs.workdir / self.subdir / f"task_{self.name}.json"
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_file, "w") as fh:
            json.dump(return_value.to_dict(), fh, indent=2)

    def result_is_valid(self, obs):
        """
        Check if the task result is valid for the observation.

        Parameters
        ----------
        obs: Observation
            The observation object.

        Returns
        -------
        bool
            True if the result is valid, False otherwise.
        """
        previous_result = self.get_result(obs)
        # currently, the cache is "valid" if it already run,
        # but we can invalidate it based on the previous return code, the time stamps of the inputs,
        # changes in the code, etc.
        valid = (
            previous_result is not None
            and previous_result.return_code != ReturnCode.INVALID
        )
        return valid

    def should_run(self, obs, requested=None, unknown_ok=True):
        """
        Check if the task should run based on a list of requested outputs.

        This function does the following checks:

            - If the stored result is invalid, the task should run.
            - If an output in the requested list does not exist, the task should run,
            - unless the task already ran succesfully and the requested output was not produced.

        Parameters
        ----------
        obs: Observation
            The observation object.
        requested: list, optional
            List of requested outputs. If None, all outputs possible are requested.
        unknown_ok: bool, optional
            If True, unknown outputs in the requested list are ignored.
            If False, an error is raised if an unknown output is requested.
            Default is True.

        Returns
        -------
        bool
            True if the task should run, False otherwise.
        """

        if not self.result_is_valid(obs):
            return True

        params = self.get_parameters(obs)

        possible_outputs = set(params["outputs"].values())

        # if not given, then request all outputs
        requested = set(possible_outputs if requested is None else requested)

        if not unknown_ok:
            # check that all the requested outputs are possible outputs of this function
            unknown_outputs = requested - possible_outputs
            if unknown_outputs:
                raise ValueError(
                    f"Unknown outputs in {self.name} task: {', '.join(unknown_outputs)}"
                )

        # ignore requests for outputs that were not produced in the previous run
        if (previous_result := self.get_result(obs)).return_code in [
            ReturnCode.OK,
            ReturnCode.SKIP,
        ]:
            requested &= set(previous_result.actual_outputs)

        return any(not obs.file_path(value).exists() for value in requested)

    def run(self, obs, requested=None, run_dependencies=False):
        """
        Run the task for the observation.

        This function does not dependency checking by default.

        Parameters
        ----------
        obs: Observation
            The observation object.
        requested: list, optional
            List of requested outputs. If None, all outputs possible are requested.
        run_dependencies: bool, optional
            If True, resolve and run dependencies and then run the task.
            If False, only run the task. It expects that the inputs exist. If an input is missing,
            an exception will be raised.

        Returns
        -------
        ReturnValue
            The return value of the task.
        """
        if run_dependencies:
            manager = self._manager()
            if manager is None:
                raise RuntimeError("Task manager is not available.")
            rv = manager.run_task(obs, self.name)
        elif self.should_run(obs, requested=requested):
            rv = self._run(obs)
        else:
            rv = self.get_result(obs)

        return rv

    def _run(self, obs):  # noqa: PLR0912
        """
        Run the task for the observation.

        This function does not dependency checking by default. It expects that the inputs exist.
        If an input is missing, an exception will be raised.
        """
        start = CxoTime.now()

        # download
        if self._download:
            obs.download(self._download)

        params = self.get_parameters(obs)

        inputs = {key: obs.file_path(val) for key, val in params["inputs"].items()}
        opt_inputs = {
            key: obs.file_path(val) for key, val in params["optional_inputs"].items()
        }
        outputs = {key: obs.file_path(val) for key, val in params["outputs"].items()}

        if set(inputs) & set(opt_inputs):
            raise ValueError(
                f"Input and optional input keys cannot overlap in task {self.name}. "
                f"Inputs: {inputs.keys()}, Optional inputs: {opt_inputs.keys()}"
            )

        for key, value in inputs.items():
            if not value.exists():
                # try with a regex if the input is a glob pattern (e.g. "primary/*_evt2.fits*")
                if match := list(obs.file_glob(params["inputs"][key])):
                    inputs[key] = match
                else:
                    inputs[key] = None

        for key, value in opt_inputs.items():
            if not value.exists():
                # try with a regex if the input is a glob pattern (e.g. "primary/*_evt2.fits*")
                if match := list(obs.file_glob(params["optional_inputs"][key])):
                    opt_inputs[key] = match
                else:
                    opt_inputs[key] = None

        if any(value is None for value in inputs.values()):
            missing_inputs = [
                str(value) for key, value in inputs.items() if value is None
            ]
            raise FileNotFoundError(
                f"Missing input files for task {self.name}: {', '.join(missing_inputs)}"
            )

        inputs.update(
            {key: value for key, value in opt_inputs.items() if value is not None}
        )

        for path in outputs.values():
            path.parent.mkdir(parents=True, exist_ok=True)

        result = self.func(obs, inputs=inputs, outputs=outputs)

        # for convenience, we allow the return value to be None,
        # a ReturnCode, or a tuple of (ReturnCode, str)
        rc = ReturnCode.OK
        msg = ""
        if result is not None:
            if isinstance(result, ReturnCode):
                rc = result
            elif isinstance(result, tuple):
                rc, msg = result
            else:
                raise TypeError(
                    f"Return value of {self.name} must be None, a ReturnCode, "
                    f"or a tuple of (ReturnCode, str), got {type(result)}"
                )

        rv = ReturnValue(
            rc,
            msg,
            start=start,
            stop=CxoTime.now(),
            inputs=params["inputs"],
            outputs=params["outputs"],
            actual_outputs=sorted(
                [
                    value
                    for value in params["outputs"].values()
                    if obs.file_path(value).exists()
                ]
            ),
        )

        self.set_result(obs, rv)

        return rv

    def skip(self, obs, msg):
        rv = ReturnValue(
            return_code=ReturnCode.SKIP,
            msg=msg,
        )
        self.set_result(obs, rv)


class TaskManager:
    """
    Class to manage tasks.
    """

    def __init__(self):
        self.tasks = {}
        self.task_graph = None

    def register_task(self, task):
        """
        Register a task.

        Parameters
        ----------
        task: Task
            The task to register.
        """
        if not isinstance(task, Task):
            raise TypeError("task must be an instance of Task")
        if task.name in self.tasks:
            raise ValueError(f"Task {task.name} is already registered.")
        self.tasks[task.name] = task
        self._set_dependency_graph()

    def run_task(self, obs, task_name):
        """
        Run a task by its name.

        Parameters
        ----------
        task_name: str
            The name of the task to run.
        obs: Observation
            The observation object.

        Returns
        -------
        ReturnValue
            A dictionary of return values for each task (requested or run) with task names as keys.
        """
        if task_name not in self.tasks:
            raise ValueError(f"Task {task_name} not found.")
        return self.run_tasks(obs, [task_name])

    def run_tasks(self, obs, task_names=None, requested_files=None):
        """
        Run task by their names.

        If any task fails, execution stops.

        Parameters
        ----------
        obs: Observation
            The observation object.
        task_names: list
            The list of names of the tasks to run.
        requested_files: list
            The list of file names to be produced by the tasks.

        Returns
        -------
        dict
            A dictionary of return values for each task (requested or run) with task names as keys.
        """
        task_names = [] if task_names is None else task_names
        requested_files = [] if requested_files is None else requested_files

        tasks = self.get_tasks_to_run(
            obs, requested_tasks=task_names, requested_files=requested_files
        )

        # requested tasks with valid stored results will not be run,
        # but their stored results will be returned because they are requested explicitly.
        return_values = {
            name: self.tasks[name].get_result(obs)
            for name in task_names
            if name not in tasks
        }

        skip = defaultdict(list)  # these will be skipped
        for name, task in tasks.items():
            for dep in self.get_dependencies(name):
                if (
                    self.tasks[dep].get_result(obs).return_code == ReturnCode.SKIP
                    and dep not in skip[name]
                ):
                    skip[name].append(dep)

            if name in skip:
                # if the task is in the skip list, it means that at least one of its dependencies
                # was skipped.
                msg = f"missing dependencies: {', '.join(skip[name])}"
                logger.info(f"Skipping task {name} for observation {obs.obsid}. {msg}")
                return_values[name] = ReturnValue(
                    return_code=ReturnCode.SKIP,
                    msg=msg,
                )
                # setting this here is a hack. Normally, the result is set inside the task,
                # but I am not running the task, and if I do not set it here, downstream tasks
                # will not be aware that this task was skipped.
                task.set_result(obs, return_values[name])
                continue

            logger.info(f"{obs} Running task {name}.")
            task.invalidate_result(obs, follow_dependents=False)
            return_value = task.run(obs)
            return_values[name] = return_value

            if return_value.return_code.value == ReturnCode.OK.value:
                logger.info(f"{obs} Task {name} completed successfully.")
            elif return_value.return_code.value == ReturnCode.SKIP.value:
                logger.info(f"{obs} Task {name} skipped: {return_value.msg}")
                for dep in self.get_dependents(name):
                    skip[dep].append(name)
            elif return_value.return_code.value == ReturnCode.WARNING.value:
                logger.warning(f"{obs} Task {name} skipped: {return_value.msg}")
                for dep in self.get_dependents(name):
                    skip[dep].append(name)
            else:
                logger.error(
                    f"{obs} Task {name} failed with return code "
                    f"{return_value.return_code.name}: {return_value.msg}"
                )
                break
        return return_values

    def _set_dependency_graph(self):
        """
        Set the dependency graph for the tasks.

        Dependencies are described by a directed acyclic graph (DAG) where nodes are tasks
        and edges represent dependencies between tasks.

        The graph is built from the inputs and outputs of each task.
        """
        task_dependency_graph = nx.DiGraph()
        for name, task in self.tasks.items():
            task_dependency_graph.add_node(name, bipartite=0)
            for fn in task._inputs.values():
                task_dependency_graph.add_node(fn, bipartite=1)
            for fn in task._optional_inputs.values():
                task_dependency_graph.add_node(fn, bipartite=1)
            for fn in task._outputs.values():
                task_dependency_graph.add_node(fn, bipartite=1)

        for name, task in self.tasks.items():
            for fn in task._inputs.values():
                task_dependency_graph.add_edge(fn, name)
            for fn in task._optional_inputs.values():
                task_dependency_graph.add_edge(fn, name)
            for fn in task._outputs.values():
                task_dependency_graph.add_edge(name, fn)

        for layer, nodes in enumerate(
            nx.topological_generations(task_dependency_graph)
        ):
            # `multipartite_layout` expects the layer as a node attribute, so add the
            # numeric layer value as a node attribute
            for node in nodes:
                task_dependency_graph.nodes[node]["layer"] = layer

        # check that this is a DAG
        if not nx.is_directed_acyclic_graph(task_dependency_graph):
            raise ValueError(
                "Task dependency graph is not a directed acyclic graph (DAG). "
                "Check the task dependencies for cycles."
            )

        # the original graph is bipartite (files are connected only tasks and tasks are connected
        # only to files). We want the task graph only.
        self.task_graph = nx.algorithms.bipartite.projected_graph(
            task_dependency_graph,
            [name for name, bp in task_dependency_graph.nodes("bipartite") if bp == 0],
        )

    def get_tasks_to_run(self, obs, requested_tasks=None, requested_files=None):
        """
        Get the tasks that need to run for the observation.

        The determination of whether a task needs to run is based on the requested tasks and files.
        The dependency tree is traversed to ensure that all dependency tasks are included.

        Parameters
        ----------
        obs: Observation
            The observation object.
        requested_tasks: list, optional
            List of task names that are requested to run.
        requested_files: list, optional
            List of file names that are requested to be produced by the tasks.

        Returns
        -------
        dict
            A dictionary of tasks that need to run, where the keys are task names and the values
            are Task objects.
        """
        unknown_tasks = set(requested_tasks) - set(self.tasks.keys())
        if unknown_tasks:
            raise ValueError(
                f"Unknown tasks requested: {', '.join(unknown_tasks)}. "
                f"Available tasks: {', '.join(self.tasks.keys())}."
            )

        requested_files = set() if requested_files is None else set(requested_files)
        requested_tasks = set() if requested_tasks is None else set(requested_tasks)

        tasks = self.task_graph

        requested_files |= {
            filename
            for name in requested_tasks
            for filename in self.tasks[name].get_parameters(obs)["outputs"].values()
        }

        # find out which tasks produce the requested files.
        # and create a dict mapping task names to filenames.
        # note that these are the filenames _after_ variable interpolation.
        params = {name: task.get_parameters(obs) for name, task in self.tasks.items()}
        task_file_map = {
            name: set(params[name]["outputs"].values()) & requested_files
            for name in self.tasks
        }
        task_file_map = {k: v for k, v in task_file_map.items() if v}

        # define the set of tasks that are possible to run.
        # This includes the requested tasks, tasks that produce the requested files,
        # and their ancestors in the dependency graph.
        possible = set(task_file_map) | requested_tasks
        possible |= {
            ancestor for task in possible for ancestor in nx.ancestors(tasks, task)
        }

        # in the following, we want to iterate over the tasks in "topological order".
        # This is the order that ensures traversal of the graph in order of dependence.
        generations = list(nx.topological_generations(tasks))

        # check if the result is valid for each possible task.
        # ---------------------------------------------------
        # Note:
        #   - we move forward in the topological order. This means that when checking each
        #     task, we have already checked its predecessors
        #     (a task needs to run if an ancestor's result is invalid).
        result_is_valid = {}
        tasks_by_generation_fwd = [
            name for gen in generations for name in gen if name in possible
        ]
        for name in tasks_by_generation_fwd:
            result_is_valid[name] = self.tasks[name].result_is_valid(obs) and all(
                result_is_valid[p] for p in tasks.predecessors(name)
            )
        # check if each task needs to run to produce the requested files.
        # ---------------------------------------------------------------
        # Note:
        #   - we move backward in the topological order. This means that first we check
        #     if a task needs to run, and if it does, we use its inputs to determine if the
        #     ancestors need to run.
        #   - requested tasks should run.
        # should_run = {task: False for task in tasks_by_generation_fwd}
        should_run = {}
        tasks_by_generation_bwd = [
            name for gen in generations[::-1] for name in gen if name in possible
        ]
        for name in tasks_by_generation_bwd:
            should_run[name] = self.tasks[name].should_run(obs, requested_files)
            if should_run[name]:
                requested_files |= set(params[name]["inputs"].values())
                requested_files |= set(params[name]["optional_inputs"].values())

        # the final result of whether each task should run is the logical OR of whether:
        #   - the stored result is invalid,
        #   - the outputs are requested.
        return {
            key: self.tasks[key]
            for key in tasks_by_generation_fwd
            if should_run[key] or (not result_is_valid[key])
        }

    def get_dependencies(self, task_name):
        """
        Get the dependencies of a task.

        Parameters
        ----------
        task_name: str
            The name of the task.
        Returns
        -------
        dict
            A dictionary of tasks on which the given task depends (the ancestors in the DAG).
        """
        return {
            name: self.tasks[name] for name in nx.ancestors(self.task_graph, task_name)
        }

    def get_dependents(self, task_name):
        """
        Get the dependents of a task.

        Parameters
        ----------
        task_name: str
            The name of the task.
        Returns
        -------
        dict
            A dictionary of tasks that depend on the given task name (the descendents in the DAG).
        """
        return {
            name: self.tasks[name]
            for name in nx.descendants(self.task_graph, task_name)
        }

    def task(
        self,
        name=None,
        inputs=None,
        optional_inputs=None,
        outputs=None,
        download=None,
        variables=None,
    ):
        """
        Decorator to run tasks.

        Parameters
        ----------
        name: str
            Name of the task.
        inputs: dict
            Dictionary of input files.
        optional_inputs: dict
            Dictionary of optional input files.
        outputs: dict
            Dictionary of output files.
        download: list
            List of file categories to download before running the task.
        variables: dict
            Dictionary of variables to be used in the task. Can be a function, with no arguments.
        """

        def wrapper(func):
            task = Task(
                func,
                manager=self,
                name=name,
                inputs=inputs,
                optional_inputs=optional_inputs,
                outputs=outputs,
                download=download,
                variables=variables,
            )
            self.register_task(task)
            return task

        return wrapper

    def dependencies(
        self,
        name=None,
        tasks=None,
        required_files=None,
        optional_files=None,
        download=None,
        variables=None,
    ):
        """
        Decorator to run tasks.

        Parameters
        ----------
        tasks: str
            Name of the task.
        required_files: list
            List of required input files. An exception is raised if any of these files are missing.
        optional_files: list
            List of optional input files.
        """

        def wrapper(func):
            task = DependentWrapper(
                func,
                manager=self,
                name=name,
                tasks=tasks,
                required_files=required_files,
                optional_files=optional_files,
                download=download,
                variables=variables,
            )
            return task

        return wrapper


class DependentWrapper(Dependent):
    """
    A wrapper for StoredResult that allows it to be used as a descriptor.
    """

    def __get__(self, obj, objtype_):
        # This method is called when the StoredResultWrapper is an attribute of an object.
        # It returns a new instance of StoredResult with the object as the "instance".
        task = Dependent(
            self.func,
            instance=obj,
            manager=self._manager,
            name=self.name,
            tasks=self._tasks,
            required_files=self._required_files,
            optional_files=self._optional_files,
            download=self._download,
            variables=self._variables,
        )
        return task


def _fullname(o):
    if hasattr(o, "name"):
        return o.name
    module = o.__class__.__module__
    if module == "builtins":
        return o.__qualname__  # avoid outputs like 'builtins.str'
    return f"{module}.{o.__qualname__}"


TASKS = TaskManager()

task = TASKS.task
dependencies = TASKS.dependencies


def run_tasks(obs, task_names=None, requested_files=None):
    """
    Run a list of tasks for an observation.

    Parameters
    ----------
    task_names: list
        List of task names to run.
    obs: Observation
        The observation object.
    """
    return TASKS.run_tasks(obs, task_names, requested_files)
