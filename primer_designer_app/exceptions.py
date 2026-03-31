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
