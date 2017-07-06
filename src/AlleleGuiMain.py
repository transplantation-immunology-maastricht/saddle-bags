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

import os

import Tkinter, Tkconstants, tkFileDialog, tkMessageBox
from Tkinter import *

from AlleleGuiEMBL import AlleleGuiEMBL
from AlleleGuiIMGT import AlleleGuiIMGT

from AlleleSubCommon import *

class AlleleGuiMain(Tkinter.Frame):

    # Initialize the GUI
    def __init__(self, root):
        Tkinter.Frame.__init__(self, root)
        root.title("An HLA Allele Submission Generator")
        self.parent = root

        self.initialize()

    # Initialize GUI elements
    def initialize(self):

        button_opt = {'fill': Tkconstants.BOTH, 'padx': 35, 'pady': 5}
        
        # Load configuration
        loadConfigurationFile()
        
        # To define the exit behavior
        self.parent.protocol('WM_DELETE_WINDOW', self.closeWindow)
        
        # This window should not be resizeable. I guess.
        self.parent.resizable(width=False, height=False)

        # Instruction Frame
        self.instructionFrame = Tkinter.Frame(self)
        self.instructionText = Tkinter.StringVar()       
        self.instructionText.set('\nSaddlebags is an HLA Allele Submission Generator.\n'
            + 'You can generate an allele submission text file for either\n'
            + 'the EMBL/ENA or IMGT/HLA nucleotide databases. You must choose:\n'
            )
        Tkinter.Label(self.instructionFrame, width=85, height=5, textvariable=self.instructionText).pack()
        self.instructionFrame.pack()
           
        # Make a frame for the more-info buttons
        self.moreInfoFrame = Tkinter.Frame(self)
        Tkinter.Button(self.moreInfoFrame, text='Generate an EMBL submission', command=lambda: self.openAlleleSubGUI('EMBL')).grid(row=0, column=0)
        Tkinter.Button(self.moreInfoFrame, text='Generate an IMGT submission', command=lambda: self.openAlleleSubGUI('IMGT')).grid(row=0, column=1)
        Tkinter.Button(self.moreInfoFrame, text='  How to use this tool   ', command=self.howToUse).grid(row=1, column=0)
        Tkinter.Button(self.moreInfoFrame, text='Contacting or Citing MUMC', command=self.contactInformation).grid(row=1, column=1)
        self.moreInfoFrame.pack()
        
        # Frame for the exit button
        self.exitFrame = Tkinter.Frame(self)
        Tkinter.Button(self.exitFrame, text='Exit', command=self.closeWindow).pack(**button_opt)
        self.exitFrame.pack()
        
        self.pack()
        
        self.initializeWindowLocation()
       
    # Put the GUI on the center of the screen. Doesn't make sense for it to start in a corner.
    # Well, lets divide by 4 instead of 2. Center is too...centered.
    def initializeWindowLocation(self):
        self.parent.update_idletasks()
        w = self.parent.winfo_screenwidth()
        h = self.parent.winfo_screenheight()
        size = tuple(int(_) for _ in self.parent.geometry().split('+')[0].split('x'))
        x = w/4 - size[0]/2
        y = h/4 - size[1]/2
        self.parent.geometry("%dx%d+%d+%d" % (size + (x, y)))
        

    # This method should popup some instruction text in a wee window.
    # This should be explicit on how to use the tool.    
    def howToUse(self):
        tkMessageBox.showinfo('How to use this tool',
            'This software is to be used to create an\n'
            + 'EMBL-formatted submission document,\n'
            + 'which specifies a (novel) HLA allele.\n\n'       
                       
            + 'This tool requires you to submit a\n'
            + 'full length HLA allele, including\n'
            + '5\' and 3\' UTRs.\n\n'
            
            + 'Use capital letters for exons,\n'
            + 'lowercase for introns & UTRs.\n\n'
            
            + 'Push the "Example Sequence" button to see a small example of'
            + ' a formatted sequence.\n'
            + 'Sequences should follow this pattern:\n'
            + '5\'utr EX1 int1 EX2 ... EX{X} 3\'utr\n\n'
            
            + 'To use this tool:\n'
            + '1.) Fill in a Sample ID, Gene Name, and Allele.'
            + ' This text will be included in the submission.\n'
            + '2.) Paste your formatted sequence in the\n'
            + 'Annotated Sequence text area.\n'
            + '3.) Push \"Generate an EMBL submission\" button'
            + ' to generate a submission.\n'
            + '4.) Push the "Save the submission" button'
            + ' to store the submission on your computer.\nYou can submit this file to EMBL.\n\n'
            
            + 'All spaces, tabs, and newlines are'
            + ' removed before the nucleotide sequence is translated.'
            )
        
        
    def contactInformation(self):
        # This method should list contact information for MUMC, and a link to the github page.  
        tkMessageBox.showinfo('Contact Information',
            'This software was created at\n'
            + 'Maastricht University Medical Center\n'
            + 'Transplantation Immunology\n'
            + 'Tissue Typing Laboratory.\n'
            + 'by Ben Matern:\n'
            + 'ben.matern@mumc.nl\n\n'
            
            + 'Please send Ben your bioinformatics\n'
            + 'and data related questions.\n\n'
            
            + 'all other inquiries can be directed\n'
            + 'to Marcel Tilanus:\n'
            + 'm.tilanus@mumc.nl\n\n'
            
            + 'This code will be hosted at:\n'
            + 'https://github.com/transplantation-\nimmunology/saddle-bags\n'
            + 'You will find more information on\n'
            + 'EMBL\'s data format on that page.'

            )
        
    def closeWindow(self):
        writeConfigurationFile()
        self.parent.destroy()        
        

    def restoreWindowPosition(self):
        # Geometry is a string that looks like this: 599x144+681+52
        # WidthxHeight+Xpos+Ypos
        newGeometry = self.windowWidth + 'x' + self.windowHeight + '+' + self.windowXpos + '+' + self.windowYpos
        self.parent.geometry(newGeometry)


    def onCloseOtherFrame(self, event):
        # <Destroy> is triggered for each widget on the subframe.
        # We want to only trigger if the main subframe is destroyed.
        if(event.widget is self.alleleSubRoot):
            self.parent.deiconify()        
            self.restoreWindowPosition()

    def rememberWindowPosition(self):
        # Remember the geometry of this window.
        self.windowWidth = str(self.parent.winfo_width())
        self.windowHeight = str(self.parent.winfo_height())
        # "Geometry" is a string that looks like this: 599x144+681+52
        # WidthxHeight+Xpos+Ypos
        windowGeometryPosTokens = self.parent.winfo_geometry().split('+')
        self.windowXpos = str(windowGeometryPosTokens[1])
        self.windowYpos = str(windowGeometryPosTokens[2])

    def openAlleleSubGUI(self, submissionType):          
        self.rememberWindowPosition()
        
        self.parent.withdraw()
        self.alleleSubRoot = Tkinter.Toplevel()
        self.alleleSubRoot.bind("<Destroy>", self.onCloseOtherFrame)
        
        if(submissionType=='IMGT'):
            print ('Opening the IMGT Submission GUI')
            AlleleGuiIMGT(self.alleleSubRoot).pack()
        elif(submissionType=='EMBL'):
            print ('Opening the EMBL Submission GUI')
            AlleleGuiEMBL(self.alleleSubRoot).pack()
        else:
            raise Exception('Unknown Submission Type.  I expected IMGT or EMBL:' + str(submissionType))
        
        # Set the X and the Y Position of the window, so it is nearby.  
        # it is necessary to update the window before assigning geometry.
        # Using Size Values from Subwindow, but Position values from Parent window 
        self.alleleSubRoot.update()
        #print('after update geometry subwindow:' + self.alleleSubRoot.winfo_geometry())
        newGeometry = (str(self.alleleSubRoot.winfo_width()) + 'x' 
            + str(self.alleleSubRoot.winfo_height()) + '+' 
            + str(self.windowXpos) + '+' 
            + str(self.windowYpos))
        self.alleleSubRoot.geometry(newGeometry)
                    
        self.alleleSubRoot.mainloop()
