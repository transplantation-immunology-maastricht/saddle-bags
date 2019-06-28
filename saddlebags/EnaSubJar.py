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

#from os.path import join

# In this file we submit to EMBL/ENA using the webin .jar file.
# ENA Submission manual can be found here:
# https://ena-docs.readthedocs.io/en/latest/general-guide/webin-cli.html

def findJarFile():
    print ('Searching for a Jar File.')
    # TODO: Put this in a config file. Hardcode the default value of just the file name, to be found inside the extracted directory.
    # Allow user to change it to a full path in the config file.
    return "/home/ben/Github/saddle-bags/jar/webin-cli-1.8.6.jar"



