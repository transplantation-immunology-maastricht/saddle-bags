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


# TODO: Check if i need all these imports. Probbably.

import sys
from sys import exc_info

try:
    from sys import _MEIPASS
    print('Running inside a compiled exe file. Great.')
except Exception:
    print('No MEIPASS Directory. This is not running from a compiled EXE file. No problem.')

from os import makedirs, name
from os.path import expanduser, join, abspath, split, isdir
from tkinter import messagebox, simpledialog

import logging

def showInfoBox(title, message):
    # A wrapper method for the tkinter popup box.
    messagebox.showinfo(title, message)

def getInfoBox(title, message):
    # wrapper to get a text input from the user.
    return simpledialog.askstring(title, message)

def showYesNoBox(title, message):
    # A wrapper method for the tkinter ask yes/no question box.
    response = messagebox.askquestion(title, message,icon='warning')
    # the response is a string, 'yes' or 'no'. That's really funny.
    return response == 'yes'

def assignIcon(tkRootWindow):
    logging.debug ('Assigning Icon for the GUI.')

    # Find window location inside executable
    try:
        # Changing this to use the join function, I don't know why I did this double slash thing before.
        # TODO: Does the icon work in the compiled exe?
        #iconFileLocation = resourcePath('images\\horse_image_icon.ico')
        windowsIconFileLocation = resourcePath(join('images', 'horse_image_icon.ico'))
        linuxIconFileLocation = resourcePath(join('images', 'horse_image_icon.xbm'))

        # Different icon code for Linux and Windows. "name"="os.name"
        if "nt" == name:
            logging.debug('I am assigning this icon:' + windowsIconFileLocation)
            tkRootWindow.wm_iconbitmap(bitmap=windowsIconFileLocation)
        else:
            logging.debug('I am assigning this icon:@' + linuxIconFileLocation)
            tkRootWindow.wm_iconbitmap(bitmap="@" + linuxIconFileLocation)
    except Exception:
        #base_path = os.path.abspath(".")
        logging.error('Could not assign icon based on path inside executable.')
        logging.error (exc_info())

    # Linux - I have given up on setting an icon in linux. I can't seem to load up any file format.
  
def resourcePath(relativePath):
    # Where will I find my resources? This should work in, or outside, a compiled EXE
    # PyInstaller creates a temp folder and stores path in _MEIPASS
    if hasattr(sys, '_MEIPASS'):
        return join(_MEIPASS, relativePath)
    return join(abspath('.'), relativePath)

def createOutputFile(outputfileName):
    # This method is a directory-safe way to open up a write file.
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        logging.info('Making Directory for output file:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput

def getSaddlebagsDirectory():
    # Do I need to detect operating system? I don't think so.
    homeDirectory = expanduser("~")
    saddlebagsDirectory = join(homeDirectory, 'saddlebags')
    return saddlebagsDirectory


