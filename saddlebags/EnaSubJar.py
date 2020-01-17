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


import logging

from os.path import join, isfile

from saddlebags.AlleleSubCommon import resourcePath
from saddlebags.SaddlebagsConfig import getConfigurationValue

# In this file we submit to EMBL/ENA using the webin .jar file.
# ENA Submission manual can be found here:
# https://ena-docs.readthedocs.io/en/latest/general-guide/webin-cli.html

def findJarFile():
    logging.debug('Searching for a Jar File.')

    configJarFileLocation = getConfigurationValue('webin_jar_location')

    if(configJarFileLocation == 'webin-cli.jar'):
        # This is the default configuration value. Normal. Use the .jar file that I have in my resource path.
        logging.debug('Using this jar file:' + str(resourcePath(join('jar', configJarFileLocation))))
        return (resourcePath(join('jar', configJarFileLocation)))
    else:
        # User has a custom jar file location in the configuration file.
        # In this case, it should be a full path to the file. Hopefully it is.
        if(isfile(configJarFileLocation)):
            logging.debug('Using this jar file:' + configJarFileLocation)
            return configJarFileLocation
        else:
            errorText = ('Invalid configuration value for jar file location. ' +
                'This is not a valid file: ' + str(configJarFileLocation))
            logging.error(errorText)
            raise Exception(errorText)





