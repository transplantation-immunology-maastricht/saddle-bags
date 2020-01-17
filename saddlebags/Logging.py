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

from os.path import join
from saddlebags.AlleleSubCommon import getSaddlebagsDirectory
import logging

def initializeLog():
    logFileLocation = join(getSaddlebagsDirectory(),'Saddlebags.Log.txt')

    # If there is no "globals" yet, then I haven't loaded a config yet, lets default to 'DEBUG' level.
    if (not ("globalVariables" in globals()) or  getConfigurationValue('logging') is None):
        logLevelText = 'DEBUG'
    else:
        logLevelText = getConfigurationValue('logging')

    logLevel = getattr(logging,logLevelText.upper())

    logFormatter = logging.Formatter("%(asctime)s:%(name)s:%(levelname)s:%(message)s")
    rootLogger = logging.getLogger()

    # Remove handlers, It's easiest for me to add my own.
    rootLogger.handlers = []
    rootLogger.setLevel(logLevel)

    # Add the Stream Handler to print to console.
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    # Add the File Handler to log to a file.
    fileHandler = logging.FileHandler(logFileLocation)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)