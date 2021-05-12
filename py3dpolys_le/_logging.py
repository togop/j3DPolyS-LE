import logging

_loggers = {}


def get_logger(name="py3dpolys_le"):
    global _loggers

    if name not in _loggers:
        _loggers[name] = logging.getLogger(name)
        _loggers[name].addHandler(logging.NullHandler())

    return _loggers[name]
