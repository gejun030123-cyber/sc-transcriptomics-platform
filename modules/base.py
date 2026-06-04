from abc import ABC, abstractmethod
from typing import Callable, Optional

class BaseAnalysis(ABC):
    MODULE_NAME = ""
    DISPLAY_NAME = ""
    DESCRIPTION = ""
    INPUT_REQUIRES = []
    PARAM_SCHEMA = {}

    def __init__(self, project_dir: str, params: dict,
                 progress_callback: Callable[[int, str], None]):
        self.project_dir = project_dir
        self.params = params
        self._progress = progress_callback

    def progress(self, pct: int, message: str):
        self._progress(min(pct, 100), message)

    @abstractmethod
    def validate_input(self, adata) -> Optional[str]:
        pass

    @abstractmethod
    def run(self, input_path: str) -> dict:
        pass
