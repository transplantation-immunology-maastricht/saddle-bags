# This file is part of saddle-bags.
#
# saddle-bags is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# saddle-bags is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with saddle-bags. If not, see <http://www.gnu.org/licenses/>.

from sys import argv, exc_info
from os.path import isfile
from os import environ
from tkinter import Tk
import subprocess
from subprocess import check_output, STDOUT, PIPE
from re import search
import logging
from saddlebags.AlleleSubCommon import showInfoBox
from saddlebags.Logging import initializeLog
from saddlebags.SaddlebagsConfig import loadConfigurationFile
from saddlebags.AlleleSubMainGui import AlleleSubMainGui
from saddlebags.EnaSubJar import findJarFile

# TODO: Version has never really been updated
SoftwareVersion = 'saddlebags Version 1.4'

def checkPrerequisites():
    logging.debug('Checking for prerequisites')

    # Do we have Java?
    # That's a complicated question. Gotta deal with lots of stuff to check that in windows, inside pyinstaller.
    try:
        # Necessary nonsense for calling command in windows.
        # Should probably move this to a common method, do that when i see a bug in embl jar file submissions.
        # True if windows:
        if hasattr(subprocess, 'STARTUPINFO'):
            logging.debug('This is Windows.')
            # On Windows, this should avoid popping up a console window, when run in --noconsole mode.
            startupInfo = subprocess.STARTUPINFO()
            startupInfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            # Pass cmd an environment so Windows will search the path variables.
            environVars = environ
            # Use an intermediate shell to launch the process? Yes, in Windows.
            useShell=True
        else:
            logging.debug('This is not Windows.')
            # we don't need these variables in linux.
            startupInfo = None
            environVars = None
            useShell = False

        # Set up some arguments for check_output
        processArgs = {
            'stdin': PIPE
            , 'stderr': STDOUT
            , 'startupinfo': startupInfo
            , 'env': environVars
            , 'universal_newlines': True
            , 'shell': useShell
        }

        txt = check_output(['java', '-version'], **processArgs)
        javaVersionOutput = str(txt)
        logging.debug('Java Version Output: ' + str(javaVersionOutput))

        # Filter out the java version. If it's there, then great.
        regexPattern = '\"(\d+\.\d+).*\"'
        javaVersion = search(regexPattern, javaVersionOutput).groups()[0]
        logging.debug('Java Version: ' + str(javaVersion))

        if (len(str(javaVersion)) < 2):
            showInfoBox('Missing Java', 'Warning.\nJava version\nwas not found.\nPerhaps java is missing?')

        # logging.debug('Java version output:\n' + javaVersionOutput)
    except Exception as e:
        showInfoBox('Missing Java', 'Warning.\nJava version\nwas not found.\nPerhaps java is missing?')
        # logging.debug ('Unexpected problem during execution:')
        logging.error('Java version was not found. Perhaps java is missing?')
        logging.debug(exc_info()[1])
        logging.debug(str(e))

    # Do i have the EMBL Commandline Jar file?
    jarFileLocation = findJarFile()
    if (isfile(jarFileLocation)):
        logging.debug('Using this EMBL Jar file:' + str(jarFileLocation))
    else:
        logging.error('This does not appear to be a valid jar file:' + str(jarFileLocation))
        showInfoBox('Missing Jar File','Warning.\nEMBL Commandline Jar File\nwas not found:\n' + str(jarFileLocation))

    # TODO: Can I see the Webservice? Somehow ping the website from Python?
    # TODO: That is, both the EMBL webservice, and the webin webservice, need em both.
    # TODO: Other Prerequisites? Should I check that important python packages are installed?
    # TODO: Anything to check for google drive submission?

if __name__=='__main__':
    try:
        # This is a really simple way to read commandline args, 
        # because there really shouldn't be any.
        # TODO: Be more graceful with this, there are better ways to read args. In fact ive written better ways.
        # No parameters are expected at all.  sys.argv[0] doesn't count.
        if (len(argv) == 1):
            print('\n\n\n\n\n\n\n\n\n\n')
            initializeLog()
            loadConfigurationFile()
            checkPrerequisites()

            logging.info('*******Starting Saddlebags*******')
            root = Tk()
            AlleleSubMainGui(root).pack()
            root.mainloop()
            logging.info('*******Closing Saddlebags*******')

            print('\n\n\n\n\n\n\n\n\n\n')

        # Print the Software Version
        elif (len(argv) == 2 and (
            argv[1].lower() == '-v' or 
            argv[1].lower() == '--version' or 
            argv[1].lower() == '-version')        
        ):
            print (SoftwareVersion)
            pass
            #

        # You executed the software wrong.  Sorry. 
        else:
            print("usage:\n" + 
                "\tRun this program using standard python call:\n" + 
                "\t$python AlleleSubmissionMain.py\n" + 
                "\tbiopython must be accessible in your python environment.  To run using Anaconda,\n"
                "\tCheck readme at https://github.com/transplantation-immunology-maastricht/saddle-bags\n"
            )


    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print ('Unexpected problem during execution:')
        print (exc_info()[1])
        raise

