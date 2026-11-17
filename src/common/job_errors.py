"""Exception types that let a failed job tell the user what to fix.

Each subclasses the builtin it replaces, so existing ``except ValueError`` /
``FileNotFoundError`` / ``CalledProcessError`` handlers keep working. Messages stay
English and may contain paths: they are for the server log, never the page.
"""
import subprocess


class ReferenceFastaError(ValueError):
    """The reference FASTA content cannot be used.

    ``reason`` is one of: invalid_characters, empty_sequence, no_sequences,
    duplicate_id, sequence_before_header, missing_id.
    """

    def __init__(self, message, reason, sequence_id=None, detail=None):
        super().__init__(message)
        self.reason, self.sequence_id, self.detail = reason, sequence_id, detail


class ReferenceMissingError(FileNotFoundError):
    """The reference FASTA no longer exists (for example removed by cleanup)."""


class ReferenceChangedError(RuntimeError):
    """The reference FASTA was replaced or edited while it was being used."""


class ToolUnavailableError(FileNotFoundError):
    """An external program (``tool``, e.g. makeblastdb or blastn) is not installed."""

    def __init__(self, tool, message):
        super().__init__(message)
        self.tool = tool


class ToolFailedError(subprocess.CalledProcessError):
    """An external program ran but did not finish its work."""

    def __init__(self, tool, returncode, cmd=None, output=None, stderr=None, message=None):
        super().__init__(returncode, cmd if cmd is not None else [tool], output, stderr)
        self.tool, self.message = tool, message

    def __str__(self):
        return self.message or super().__str__()
