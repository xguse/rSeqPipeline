import subprocess
import os
import sys
from rSeq.utils.errors import *

# ++++++++ Verifiying/preparing external environment ++++++++
def whereis(program):
    for path in os.environ.get('PATH', '').split(':'):
        if os.path.exists(os.path.join(path, program)) and not os.path.isdir(os.path.join(path, program)):
            return os.path.join(path, program)
    return None

def mkdirp(path):
    """Create new dir while creating any parent dirs in the path as needed.
    """

    if not os.path.isdir(path):
        try:
            os.makedirs(path)
        except OSError as errTxt:
            if "File exists" in errTxt:
                sys.stderr.write("FYI: %s" % (errTxt))
            else:
                raise
            
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def runExternalApp(progName,argStr):
    """Convenience func to handle calling and monitoring output of external programs."""
    
    # Ensure program is callable.
    if not whereis(progName):
        raise SystemCallError(None,'"%s" command not found in your PATH environmental variable.' % (progName))
    
    # Construct shell command
    cmdStr = "%s %s" % (progName,argStr)
    
    # Set up process obj
    process = subprocess.Popen(cmdStr,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    # Get results
    result  = process.communicate()
    
    # Check returncode for success/failure
    if process.returncode != 0:
        raise SystemCallError(process.returncode,result[1],progName)
    
    # Return result
    return result
