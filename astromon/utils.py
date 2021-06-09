import os
import functools
import logging
from pathlib import Path
from contextlib import AbstractContextManager


def set_ciao_context(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    os.environ["ASCDS_WORK_PATH"] = str(directory)
    pf = "{};{}:{}".format(
        str(directory.absolute()),
        os.environ["ASCDS_INSTALL"] + "/param",
        os.environ["ASCDS_INSTALL"] + "/contrib/param"
    )
    os.environ["PFILES"] = pf


class chdir(AbstractContextManager):
    def __init__(self, directory):
        self.cwd = os.getcwd()
        self.directory = directory
        os.chdir(self.directory)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        os.chdir(self.cwd)
        return exc_type is None


class CIAOContext(AbstractContextManager):
    def __init__(self, directory, logger=logging):
        self.directory = directory
        self.logger = logger
        self.logger.debug(f'entering ciao context {self.directory}')
        self.env = {k: os.environ.get(k, None) for k in ["ASCDS_WORK_PATH", "PFILES"]}
        set_ciao_context(self.directory)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        os.environ.update({
            key: val for key, val in self.env.items() if val is not None
        })
        self.logger.debug(f'leaving ciao context {self.directory}')
        return exc_type is None


class CIAOContextDecorator:
    def __init__(self, directory, logger=logging):
        self.directory = directory
        self.logger = logger

    def __call__(self, func):
        @functools.wraps(func)
        def ciao_context_wrapper(*args, **kwargs):
            with CIAOContext(self.directory, logger=self.logger):
                return func(*args, **kwargs)
        return ciao_context_wrapper


def ciao_context_function(func):
    """
    Decorator to set the CIAO parameters directory within the context of a member function.

    By default, CIAO stores function parameters in $HOME/cxcds_param4. Within a function decorated
    with this decorator, this changed to ${self.workdir}/${self.obsid}/params by setting the PFILES
    variable. This also sets ASCDS_WORK_PATH.

    The requirement is that the decorated class must have a 'workdir' and am 'obsid' data members.
    """
    @functools.wraps(func)
    def ciao_context_wrapper(self, *args, **kwargs):
        with CIAOContext(
                self.workdir / self.obsid / 'params',
                logger=logging.getLogger('astromon.utils')):
            return func(self, *args, **kwargs)
    return ciao_context_wrapper


ciao_context = CIAOContext


class LoggingCallDecorator:
    def __init__(
        self,
        logger=None,
        log_level=logging.DEBUG,
        ignore_exceptions=False
    ):
        self.logger = logger if logger is not None else logging
        self.log_level = log_level
        self.ignore_exceptions = ignore_exceptions

    def __call__(self, func):
        @functools.wraps(func)
        def _logging_call_decorator(*args, **kwargs):
            args_str = []
            if args:
                args_str += [', '.join([str(arg) for arg in args])]
            if kwargs:
                args_str += [', '.join([f'{k}={v}' for k, v in kwargs.items()])]
            args_str = ', '.join(args_str)
            self.logger.log(
                self.log_level,
                f"{func.__name__}({args_str}) Started"
            )
            try:
                result = func(*args, **kwargs)
                self.logger.log(
                    self.log_level,
                    f"{func.__name__}({args_str}) Finished"
                )
                return result
            except Exception as e:
                self.logger.log(
                    self.log_level,
                    f"{func.__name__}({args_str}) Error: {e}"
                )
                if not self.ignore_exceptions:
                    raise
        return _logging_call_decorator


logging_call_decorator = LoggingCallDecorator(
    log_level=logging.INFO,
    logger=logging.getLogger('astromon')
)
