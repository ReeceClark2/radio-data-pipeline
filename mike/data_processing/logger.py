# Standard library
import io
import logging
from typing import Dict, List


class Log_Collector:
    def __init__(self, name: str = "radio_pipeline"):
        self.log_stream = io.StringIO()
        self.logger = logging.getLogger(name)
        self.logger.setLevel(logging.DEBUG)
        self.log_handler = logging.StreamHandler(self.log_stream)
        self.log_handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        if not any(isinstance(h, logging.StreamHandler) and h.stream == self.log_stream
           for h in self.logger.handlers):
                self.logger.addHandler(self.log_handler)
        self.logger.propagate = False

    def get_log_text(self) -> str:
        return self.log_stream.getvalue()

    def get_log_entries(self) -> List[Dict[str, str]]:
        lines = self.get_log_text().strip().splitlines()
        entries = []
        for line in lines:
            if line.startswith("["):
                level, message = line.split("] ", 1)
                entries.append({
                    "level": level.strip("[]"),
                    "message": message
                })
        return entries

    # Forward logging methods to internal logger
    def debug(self, msg: str, *args, **kwargs):
        self.logger.debug(msg, *args, **kwargs)

    def info(self, msg: str, *args, **kwargs):
        self.logger.info(msg, *args, **kwargs)

    def warning(self, msg: str, *args, **kwargs):
        self.logger.warning(msg, *args, **kwargs)

    def error(self, msg: str, *args, **kwargs):
        self.logger.error(msg, *args, **kwargs)

    def save(self, filepath: str = "log_output.txt"):
        if filepath.endswith(".fits"):
            filepath = filepath[:-5] + "_log.txt"  # remove last 5 chars ('.fits') and add '.txt'
        with open(filepath, "w") as f:
            f.write(self.get_log_text())