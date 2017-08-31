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

import Tkinter
import sys

# This exception handling is necessary for packaging the modules in pyinstaller.
# The namespace(?) is required to find the modules in pyinstaller.  
try:
    from AlleleSubMainGui import AlleleSubMainGui
    from AlleleSubCommon import loadConfigurationFile
except:
    from saddlebags.AlleleSubMainGui import AlleleSubMainGui
    from saddlebags.AlleleSubCommon import loadConfigurationFile
    
SoftwareVersion = 'saddlebags Version 1.1'
    
if __name__=='__main__':
    try:
        # This is a really simple way to read commandline args, 
        # because there really shouldn't be any.
        # TODO: Be more graceful with this, there are better ways to read args.
        # No parameters are expected at all.  sys.argv[0] doesn't count.
        if (len(sys.argv) == 1):
            
            loadConfigurationFile()
            
            print('\n\n\n\n\n***Starting the HLA Allele Submission Tool***\n')

            root = Tkinter.Tk()
            AlleleSubMainGui(root).pack()
            root.mainloop()

            print('Done.  Hooray.')

        # Print the Software Version
        elif (len(sys.argv) == 2 and (
            sys.argv[1].lower() == '-v' or 
            sys.argv[1].lower() == '--version' or 
            sys.argv[1].lower() == '-version')        
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
                "\tCheck readme at https://github.com/transplantation-immunology/saddle-bags\n"
            )


    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print 'Unexpected problem during execution:'
        print sys.exc_info()[1]
        raise