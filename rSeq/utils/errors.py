
class UnexpectedValueError(Exception):
    """When values that "should" not be possible happen; like if a variable was changed unexpectedly."""
    pass

class InvalidFileFormatError(Exception):
    """When errors occur due to malformed file formats."""
    pass

class InvalidOptionError(Exception):
    pass

class ExternalError(EnvironmentError):
    """Exception raised when a problem occurs while attempting to run an external system call."""
    def __str__(self):
        if not self.filename: 
            return """ERROR: %s.\nRETURN_STATE: %s.""" % (self.strerror.strip('\n'),
                                                          self.errno)
        else: 
            return """ERROR in %s: %s.\nRETURN_STATE: %s.""" % (self.filename,
                                                                self.strerror.strip('\n'),
                                                                self.errno)

    
    

    