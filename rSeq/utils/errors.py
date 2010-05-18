
class UnexpectedValueError(Exception):
    pass

class InvalidFileFormatError(Exception):
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

    
    

    