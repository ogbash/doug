#!/usr/bin/env python

"""Common classes used by all DOUG scripts"""

class ScriptException (Exception):
    """Script running error that is understood by script, ie failure running some programs.

    Some day later, additional information can be added to this class that makes possible logs,
    error files, input files or other information that can help understand error source."""

    def __init__(self, descr):
        Exception.__init__(self, descr)
        self.files = []

    def addFile(self, fname, description=None):
        self.files.append((fname, description))
