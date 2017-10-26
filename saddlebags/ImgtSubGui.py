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

#import os
from os.path import expanduser

import Tkinter, Tkconstants, tkFileDialog, tkMessageBox
from Tkinter import Scrollbar, BOTTOM, RIGHT, X, Y, HORIZONTAL, NONE, NORMAL, DISABLED

from ImgtSubGenerator import ImgtSubGenerator
from ImgtSubOptionsForm import ImgtSubOptionsForm
from AlleleSubCommon import collectAndValidateRoughSequence, assignIcon, assignConfigurationValue, parseExons, getConfigurationValue, identifyGenomicFeatures, isSequenceAlreadyAnnotated, collectRoughSequence, getAlleleDescription

from AlleleSubCommonRest import fetchSequenceAlleleCallWithGFE

from HlaSequenceException import HlaSequenceException


#from HLAGene import HLAGene

# The AlleleGui class is an extension of Tkinter.  The GUI elements and interactions are specified in this class.
class ImgtSubGui(Tkinter.Frame):

    # I shouldn't need to write a select-All method but TK is kind of annoying.
    def selectall(self, event):

        event.widget.tag_add("sel","1.0","end")
        
        
    # Initialize the GUI
    def __init__(self, root):
        Tkinter.Frame.__init__(self, root)
        root.title("Create and Save an IPD-IMGT/HLA Sequence Submission")
        self.parent = root

        # Assign the icon of this sub-window.        
        assignIcon(self.parent)

        # Ctrl-A doesn't work by default in TK.  I guess I need to do it myself.
        root.bind_class("Text","<Control-a>", self.selectall)
        
        # To define the exit behavior.  Save the input sequence text.
        self.parent.protocol('WM_DELETE_WINDOW', self.saveAndExit)

        button_opt = {'fill': Tkconstants.BOTH, 'padx': 35, 'pady': 5}
        
        # A frame for the Instructions Label.
        self.instructionsFrame = Tkinter.Frame(self)  
        self.instructionText = Tkinter.StringVar()       
        self.instructionText.set('\nThis tool will generate an HLA allele submission for\n'
            + 'the IPD-IMGT/HLA nucleotide database.\n'
            + 'For more information:\n')
        Tkinter.Label(self.instructionsFrame, width=85, height=6, textvariable=self.instructionText).pack()
        self.instructionsFrame.pack(expand=False, fill='both')
        
        # Make a frame for the more-info buttons
        self.moreInfoFrame = Tkinter.Frame(self)
        self.howToUseButton = Tkinter.Button(self.moreInfoFrame, text='How to use this tool', command=self.howToUse)
        self.howToUseButton.grid(row=0, column=0)
        self.exampleButton = Tkinter.Button(self.moreInfoFrame, text='Example Sequence', command=self.sampleSequence)
        self.exampleButton.grid(row=0, column=1)
        self.moreInfoFrame.pack() 
       
        # Create a frame for the input widget, add scrollbars.
        self.featureInputFrame = Tkinter.Frame(self)
        
        self.featureInstrText = Tkinter.StringVar()
        self.featureInstrText.set('Annotated Sequence:')
        self.featureInstrLabel = Tkinter.Label(self.featureInputFrame, width=80, height=1, textvariable=self.featureInstrText).pack()

        self.featureInputXScrollbar = Scrollbar(self.featureInputFrame, orient=HORIZONTAL)
        self.featureInputXScrollbar.pack(side=BOTTOM, fill=X)

        self.featureInputYScrollbar = Scrollbar(self.featureInputFrame)
        self.featureInputYScrollbar.pack(side=RIGHT, fill=Y)

        self.featureInputGuiObject = Tkinter.Text(
            self.featureInputFrame
            , width=80, height=8
            , wrap=NONE
            , xscrollcommand=self.featureInputXScrollbar.set
            , yscrollcommand=self.featureInputYScrollbar.set
        )

        self.featureInputXScrollbar.config(command=self.featureInputGuiObject.xview)
        self.featureInputYScrollbar.config(command=self.featureInputGuiObject.yview) 

        self.featureInputGuiObject.pack(expand=True, fill='both') 
        self.featureInputFrame.pack(expand=True, fill='both')


        # Create  Frame for "Generate Submission" button.
        self.submButtonFrame = Tkinter.Frame(self)
        self.submissionOptionsButton = Tkinter.Button(self.submButtonFrame, text='1) Submission Options', command=self.chooseSubmissionOptions)
        self.submissionOptionsButton.grid(row=0, column=0)
        self.annotateFeaturesButton = Tkinter.Button(self.submButtonFrame, text='2) Annotate Exons & Introns' , command=self.annotateInputSequence)
        self.annotateFeaturesButton.grid(row=0, column=1)
        self.generateSubmissionButton = Tkinter.Button(self.submButtonFrame, text='3) Generate an IMGT/HLA submission', command=self.constructSubmission)
        self.generateSubmissionButton.grid(row=0, column=2)
        self.submButtonFrame.pack()
       
        # Output interface is contained on a frame.
        self.submOutputFrame = Tkinter.Frame(self)
        
        self.outputEMBLSubmission = Tkinter.StringVar()
        self.outputEMBLSubmission.set('Allele Submission Preview:')
        self.outputEMBLLabel = Tkinter.Label(self.submOutputFrame, width=80, height=1, textvariable=self.outputEMBLSubmission).pack()

        self.submOutputXScrollbar = Scrollbar(self.submOutputFrame, orient=HORIZONTAL)
        self.submOutputXScrollbar.pack(side=BOTTOM, fill=X)

        self.submOutputYScrollbar = Scrollbar(self.submOutputFrame)
        self.submOutputYScrollbar.pack(side=RIGHT, fill=Y)

        self.submOutputGuiObject = Tkinter.Text(
            self.submOutputFrame, width=80, height=8, wrap=NONE
            , xscrollcommand=self.submOutputXScrollbar.set
            , yscrollcommand=self.submOutputYScrollbar.set
        )

        self.submOutputXScrollbar.config(command=self.submOutputGuiObject.xview)
        self.submOutputYScrollbar.config(command=self.submOutputGuiObject.yview) 

        self.submOutputGuiObject.pack(expand=True, fill='both') 
        self.submOutputFrame.pack(expand=True, fill='both')

        self.uploadSubmissionFrame = Tkinter.Frame(self)        
        self.saveSubmissionButton = Tkinter.Button(self.uploadSubmissionFrame, text='4) Save Submission to My Computer', command=self.saveSubmissionFile)
        self.saveSubmissionButton.pack(**button_opt)
        # TODO: Enable IMGT Upload functionality. Also enable this button.
        self.uploadButton = Tkinter.Button(self.uploadSubmissionFrame, text='5) Upload Submission to IPD-IMGT/HLA', command=self.uploadSubmission, state=DISABLED)        
        self.uploadButton.pack(**button_opt)
        self.exitButton = Tkinter.Button(self.uploadSubmissionFrame, text='Exit', command=self.saveAndExit)
        self.exitButton.pack(**button_opt)
        self.uploadSubmissionFrame.pack()
        
        self.pack(expand=True, fill='both')

    def chooseSubmissionOptions(self):
        print ('Opening the IMGT/HLA Submission Options Dialog')
        
        self.disableGUI()
        
        imgtOptionsRoot = Tkinter.Toplevel()
        imgtOptionsRoot.bind("<Destroy>", self.enableGUI)
        ImgtSubOptionsForm(imgtOptionsRoot).pack()

        # Set the X and the Y Position of the options window, so it is nearby.  
        imgtOptionsRoot.update()        
        windowXpos = str(self.parent.winfo_geometry().split('+')[1])
        windowYpos = str(self.parent.winfo_geometry().split('+')[2])
        newGeometry = (str(imgtOptionsRoot.winfo_width()) + 'x' 
            + str(imgtOptionsRoot.winfo_height()) + '+' 
            + str(windowXpos) + '+' 
            + str(windowYpos))
        imgtOptionsRoot.geometry(newGeometry)
        
        imgtOptionsRoot.mainloop()

        
    def sampleSequence(self):
        self.featureInputGuiObject.delete('1.0','end')
        self.featureInputGuiObject.insert('1.0', 'aag\nCGTCGT\nccg\nGGCTGA\naat')
        
        # Clear the password, keep the username
        assignConfigurationValue('imgt_password','')
        
        assignConfigurationValue("allele_name",'Allele:01:02')
        assignConfigurationValue('gene','HLA-C')
        assignConfigurationValue('sample_id', 'Donor_12345')
        assignConfigurationValue('class','1')
        
        assignConfigurationValue('embl_sequence_accession', 'LT123456')        
        assignConfigurationValue('embl_release_date', '01/01/2020')
        
        assignConfigurationValue('is_published','0')
        
        assignConfigurationValue('reference_title', 'Published Reference Title') 
        assignConfigurationValue('reference_authors', 'Albert Authorman, Ben Bioinformaticist, Cindy Cell-Culture') 
        assignConfigurationValue('reference_journal', 'Scientific Journal of Research')  
        
        assignConfigurationValue('closest_known_allele', 'HLA-C*01:02:01')
        assignConfigurationValue('closest_allele_written_description', 'This allele has a C->G polymorphism in Exon 1.\nPosition 5 in the coding sequence.\nThis polymorphism is interesting because of science.')
        
        assignConfigurationValue('ethnic_origin', 'Unknown')
        assignConfigurationValue('sex', 'Unknown')
        assignConfigurationValue('consanguineous', 'Unknown')
        assignConfigurationValue('homozygous', 'Unknown')

        
        self.constructSubmission()
        
    # This method should popup some instruction text in a wee window.
    # This should be explicit on how to use the tool.    
    def howToUse(self):
        tkMessageBox.showinfo('How to use this tool',
            'This software is to be used to create an\n'
            + 'IPD-IMGT/HLA-formatted submission document,\n'
            + 'which specifies a (novel) HLA allele.\n\n'       
                       
            + 'This tool requires you to submit a\n'
            + 'full length HLA allele, including\n'
            + '5\' and 3\' UTRs.\n\n'
     
            + 'To create an IPD-IMGT/HLA submission:\n\n'
            + '1.) Paste a full-length HLA sequence in\n'
            + 'the Annotated Sequence text area.\n'
            + '2.) Push [Submission Options] and provide\n'
            + 'the necessary sequence metadata.\n'
            + '3.) Push [Annotate Exons & Introns] to\n'
            + 'annotate your exons automatically.\n'
            + '4.) Push [Generate an IPD-IMGT/HLA submission]\n'
            + 'button to generate a submission.\n'
            # TODO: Change these instructions when the submission functionality works.
            #+ '5.) Push [Upload Submission to EMBL]\n'
            #+ 'to submit the sequence\n'
            #+ 'using EMBL Webin REST interface\n\n'
            
            + 'If exon annotation is not available,\n'
            + 'it may be necessary to annotate manually.\n\n'

            + 'Sequences should follow this pattern:\n'
            + '5\'utr EX1 int1 EX2 ... EX{X} 3\'utr\n\n'  
            
            + 'Use capital letters for exons,\n'
            + 'lowercase for introns & UTRs.\n\n'
            
            + 'Push the "Example Sequence" button to see\n'
            + 'an example of a formatted sequence.\n\n'
            
            + 'More information available\n'
            + 'on the MUMC Github Page:\n'
            + 'https://github.com/transplantation-\n'
            + 'immunology-maastricht/saddle-bags'

            )
        
    # Ask user for a output file location, and write the IMGT submission to a file.
    # This takes the input from the output field, rather than generate a new submission.
    # So the user can edit the submission before or after saving it.
    def saveSubmissionFile(self):

        self.dir_opt = options = {}
       
        options['initialdir'] = expanduser("~")
        options['parent'] = self
        options['title'] = 'Specify your output file.'
        options['initialfile'] = 'IMGT.HLA.Submission.txt'
        outputFileObject = tkFileDialog.asksaveasfile(**self.dir_opt)
        submissionText = self.submOutputGuiObject.get('1.0', 'end')
        outputFileObject.write(submissionText)
        # TODO: Did I detect any exceptions? Maybe I don't have permission to write that file        
        # I saw an error when i wrote to a network drive once. 
        
    def annotateInputSequence(self): 
        try:
            self.disableGUI()
            self.update()   
            
            # Popup.  This uses NMDP BTM ACT tool to annotate sequences.
            if (tkMessageBox.askyesno('Annotate Exons?'
                , 'This will annotate your exons using the\n'
                + 'NMDP: BeTheMatch Gene Feature\n'
                + 'Enumeration / Allele Calling Tool.\n\n'
                + 'Please verify you have chosen the correct\n'
                + 'HLA Gene Locus in the\n' 
                + 'Submission Options menu.\n\n'
                + 'Do you want to continue?'                
                )):

                roughNucleotideSequence = collectAndValidateRoughSequence(self.featureInputGuiObject)
            
                alleleCallWithGFE = fetchSequenceAlleleCallWithGFE(roughNucleotideSequence, getConfigurationValue('gene'))
                annotatedSequence = parseExons(roughNucleotideSequence, alleleCallWithGFE)
                self.featureInputGuiObject.delete('1.0','end')    
                self.featureInputGuiObject.insert('1.0', annotatedSequence) 
                
                # TODO: I have not implemented this yet.    
                alleleDescription = getAlleleDescription()
                fdfsgdf                
    
            self.update()
            self.enableGUI()
            
        except Exception, e:
            tkMessageBox.showerror('Error Annotating Input Sequence.'
                , str(e))
            raise
        
    def uploadSubmission(self):
        # TODO: Implement this method.
        print ('This functionality is disabled until it works better.')


    # Gather sequence information from the input elements, and generate a text IMGT submission.
    #def constructSubmission(self):
    #    try:
    #        allGen = ImgtSubGenerator()
    #        roughFeatureSequence = self.featureInputGuiObject.get('1.0', 'end')

    #        allGen.sequenceAnnotation = identifyGenomicFeatures(roughFeatureSequence)
    #        imgtSubmission = allGen.buildIMGTSubmission()
            
    #        if (imgtSubmission is None or len(imgtSubmission) < 1):
    #            tkMessageBox.showerror('Empty submission text'
    #                ,'You are missing some required information.\n'
    #                + 'Try the \'Submission Options\' button.\n')
                
    #            self.submOutputGuiObject.delete('1.0','end')    
    #            self.submOutputGuiObject.insert('1.0', '') 
    #        else:
    #            self.submOutputGuiObject.delete('1.0','end')    
    #            self.submOutputGuiObject.insert('1.0', imgtSubmission) 
            
            
    #    except KeyError, e:
    #        tkMessageBox.showerror('Missing Submission Options'
    #            ,'You are missing some required information.\n'
    #            + 'Use the \'Submission Options\' button.\n'
    #            + 'Missing Data: ' + str(e))
                

            
    # Gather sequence information from the input elements, and generate a text EMBL submission.
    def constructSubmission(self):
        try:
        
            roughNucleotideSequence = collectAndValidateRoughSequence(self.featureInputGuiObject)
            
            if (isSequenceAlreadyAnnotated(roughNucleotideSequence)):
                annotatedSequence = roughNucleotideSequence
                
            else:
                if (tkMessageBox.askyesno('Auto - Annotate Exons?'
                    , 'It looks like your sequence features have not been identified.\n' +
                    'Would you like to annotate using NMDP: BeTheMatch\n' +
                    'Gene Feature Enumeration Tool?')):
                    
                    self.annotateInputSequence()
                    annotatedSequence = collectRoughSequence(self.featureInputGuiObject)
                else:
                    # You chose not to annotate.  Hope this works out for you.
                    annotatedSequence = roughNucleotideSequence
                
            allGen = ImgtSubGenerator()
            allGen.sequenceAnnotation = identifyGenomicFeatures(annotatedSequence)
            
            imgtSubmission = allGen.buildIMGTSubmission()
                        
            if (imgtSubmission is None or len(imgtSubmission) < 1):
                tkMessageBox.showerror('Empty submission text'
                    ,'You are missing some required information.\n'
                    + 'Try the \'Submission Options\' button.\n')
                
                self.submOutputGuiObject.delete('1.0','end')    
                self.submOutputGuiObject.insert('1.0', '') 
            else:
                self.submOutputGuiObject.delete('1.0','end')    
                self.submOutputGuiObject.insert('1.0', imgtSubmission) 
            
        except KeyError, e:
            tkMessageBox.showerror('Missing Submission Options'
                ,'You are missing some required information.\n'
                + 'Use the \'Submission Options\' button.\n'
                + 'Missing Data: ' + str(e))
           
        except HlaSequenceException, e:
            tkMessageBox.showerror('I see a problem with Sequence Format.'
                , str(e))
           
        except Exception, e:
            tkMessageBox.showerror('Error Constructing Submission.'
                , str(e))
            raise       
            
            
            
    def saveAndExit(self):
        assignConfigurationValue('sequence', self.featureInputGuiObject.get('1.0', 'end'))
        self.parent.destroy()
        
        
    def enableGUI(self, event=None):
        self.toggleGUI(True)  
        
    def disableGUI(self):
        self.toggleGUI(False)   
        
    def toggleGUI(self, isEnabled): 
        #print ('Toggling GUI Widgets:' + str(isEnabled))
         
        newState = (NORMAL if (isEnabled) else DISABLED)
        
        # Choosing the widgets individually, this makes the most sense I think.
        self.howToUseButton.config(state=newState) 
        self.exampleButton.config(state=newState)         
        #self.featureInputGuiObject.config(state=newState)
        self.submissionOptionsButton.config(state=newState)
        self.generateSubmissionButton.config(state=newState)
        self.annotateFeaturesButton.config(state=newState)
        self.submOutputGuiObject.config(state=newState)
        #self.uploadButton.config(state=newState)
        self.saveSubmissionButton.config(state=newState)
        self.exitButton.config(state=newState)
        
            

