class PrimerDesignerError(Exception):
    pass


class InvalidTranscriptIdError(PrimerDesignerError):
    pass


class InvalidTranscriptVersionError(PrimerDesignerError):
    pass


class InvalidTranscriptInputError(PrimerDesignerError):
    pass


class ExonExonJunctionError(PrimerDesignerError):
    pass


class InvalidReferenceSequenceError(PrimerDesignerError):
    """Reference template missing or unsuitable (e.g. Ensembl soft-masked to N)."""


class NoPrimerPairsFoundError(PrimerDesignerError):
    """Primer3 completed but returned zero primer pairs."""
