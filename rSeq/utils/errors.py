class SanityCheckError(StandardError):
    """When values that "should" not be possible happen; like if a variable was changed unexpectedly."""
    pass

class UnexpectedValueError(StandardError):
    """When values that "should" not be possible happen; like if a variable was changed unexpectedly."""
    pass

class InvalidFileFormatError(StandardError):
    """When errors occur due to malformed file formats."""
    pass

class MissingArgumentError(StandardError):
    """When a required argument is missing from the parsed command line options."""
    def __init__(self,errMsg):
        self.msg = errMsg
    def __str__(self):
        return """ERROR: %s""" % (self.msg)

class InvalidOptionError(StandardError):
    def __init__(self,optVal,optName,validVals=None):
        self.optVal    = optVal
        self.optName   = optName
        self.validVals = validVals
        
    def __str__(self):
        if self.validVals:
            return """ERROR: %s is not a valid value for arg:%s.\n\tValid values are: %s""" % (self.optVal,self.optName,self.validVals)
        else:
            return """ERROR: %s is not a valid value for arg:%s.""" % (self.optVal,self.optName)

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

    
    

    