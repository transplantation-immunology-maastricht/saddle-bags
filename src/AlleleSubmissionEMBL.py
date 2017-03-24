# This file is part of EMBL-HLA-Submission.
#
# EMBL-HLA-Submission is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# EMBL-HLA-Submission is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with EMBL-HLA-Submission. If not, see <http://www.gnu.org/licenses/>.

# Version 1.0 

SoftwareVersion = "EMBL-HLA-Submission Version 1.0"

import Tkinter
import sys

from AlleleGui import AlleleGui

if __name__=='__main__':
    try:
        # This is a really simple way to read commandline args, 
        # because there really shouldn't be any.
        # TODO: Be more graceful with this, there are better ways to read args.

        # No parameters are expected at all.  sys.argv[0] doesn't count.
        if (len(sys.argv) == 1):
            print('\n\n\n\n\n***Creating an EMBL Allele submission***\n')

            root = Tkinter.Tk()
            AlleleGui(root).pack()
            root.mainloop()

            print('Done.  Yay.')

        # Print the Software Version
        elif (len(sys.argv) == 2 and (
            sys.argv[1].lower() == '-v' or 
            sys.argv[1].lower() == '--version' or 
            sys.argv[1].lower() == '-version')        
        ):
            print (SoftwareVersion)

        # You executed the software wrong.  Sorry. 
        else:
            print("usage:\n" + 
                "\tRun this program using standard python call:\n" + 
                "\t$python AlleleSubmissionEMBL.py\n" + 
                "\tbiopython must be accessible in your python environment.  To run using Anaconda,\n"
                "\tCheck readme at https://github.com/transplantation-immunology/EMBL-HLA-Submission\n"
            )


    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print 'Unexpected problem during execution:'
        print sys.exc_info()[1]
        raise

