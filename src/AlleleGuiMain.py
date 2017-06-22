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

SoftwareVersion = "Bhast Version 1.1"

import os

import Tkinter, Tkconstants, tkFileDialog, tkMessageBox
from Tkinter import *

from SubmissionGeneratorEMBL import SubmissionGeneratorEMBL
from HLAGene import *

from AlleleGuiEMBL import AlleleGuiEMBL
from AlleleGuiIMGT import AlleleGuiIMGT

from AlleleSubCommon import *


# TODO: I should hide this window when the user selects an option.
# Example maybe like this:
# https://www.blog.pythonlibrary.org/2012/07/26/tkinter-how-to-show-hide-a-window/


# The AlleleGui class is an extension of Tkinter.  The GUI elements and interactions are specified in this class.
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

        # Instruction Frame
        self.instructionFrame = Tkinter.Frame(self)
  
        self.instructionText = Tkinter.StringVar()       
        self.instructionText.set('\nSaddlebags is an HLA Allele Submission Generator.\n'
            + 'You can generate an allele submission text file for either\n'
            + 'the EMBL/ENA or IMGT/HLA nucleotide databases. You must choose:\n'
            )
        Tkinter.Label(self, width=85, height=5, textvariable=self.instructionText).pack()
           
        # Submit Button Frame
        # EMBL Submission Button 
        self.submitButtonFrame = Tkinter.Frame(self)  
        Tkinter.Button(self, text='Generate an EMBL submission', command=self.openEMBLGUI).pack(**button_opt)
        Tkinter.Button(self, text='Generate an IMGT submission', command=self.openIMGTGUI).pack(**button_opt)
        
        # Make a frame for the more-info buttons
        self.moreInfoFrame = Tkinter.Frame(self)
  
        Tkinter.Button(self.moreInfoFrame, text='  How to use this tool   ', command=self.howToUse).grid(row=0, column=0)
        Tkinter.Button(self.moreInfoFrame, text='Contacting or Citing MUMC', command=self.contactInformation).grid(row=0, column=1)
        
        self.moreInfoFrame.pack()
        
        # Frame for the exit button
        self.exitFrame = Tkinter.Frame(self)
        Tkinter.Button(self.exitFrame, text='Exit', command=self.closeWindow).pack(**button_opt)
        self.exitFrame.pack()

        

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
            + 'https://github.com/transplantation-\nimmunology/EMBL-HLA-Submission\n'
            + 'You will find more information on\n'
            + 'EMBL\'s data format on that page.'

            )
        
        
        
        
    def closeWindow(self):
        writeConfigurationFile()

        self.parent.destroy()        
        

        
    # TODO: this isn't workng right now, i have to set handlers and such.
    # https://www.blog.pythonlibrary.org/2012/07/26/tkinter-how-to-show-hide-a-window/
    # A method to call when the subframe is closed.
    def onCloseOtherFrame(self, otherFrame):
        otherFrame.destroy()
        #elf.show() 
        self.parent.update()
        self.parent.deiconify()   

    def openEMBLGUI(self):
        print ('Opening the EMBL Submission GUI')
        # TODO: Uncomment this when the subwindows are all working right.
        #self.parent.withdraw()
        emblSubRoot = Tkinter.Toplevel()
        AlleleGuiEMBL(emblSubRoot).pack()
        emblSubRoot.mainloop()
        
    def openIMGTGUI(self):
        print ('Opening the IMGT Submission GUI')
        imgtSubRoot = Tkinter.Toplevel()
        AlleleGuiIMGT(imgtSubRoot).pack()
        imgtSubRoot.mainloop()
