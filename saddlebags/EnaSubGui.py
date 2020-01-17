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

from sys import exc_info

from os.path import expanduser

from tkinter import Radiobutton, Entry, Menu, Frame, StringVar, IntVar, Button, Text, Label, Toplevel, Scrollbar, messagebox, filedialog
from tkinter.constants import BOTH, BOTTOM, RIGHT, X, Y, NONE, HORIZONTAL, NORMAL, DISABLED

from saddlebags.EnaSubGenerator import EnaSubGenerator
from saddlebags.EnaSubOptionsForm import EnaSubOptionsForm
from saddlebags.AlleleSubCommon import assignIcon, showInfoBox
from saddlebags.HlaSequence import collectAndValidateRoughSequence
from saddlebags.SaddlebagsConfig import writeConfigurationFile, getConfigurationValue, assignConfigurationValue
from saddlebags.AlleleSubmission import AlleleSubmission, SubmissionBatch
from saddlebags.HlaSequenceException import HlaSequenceException

import logging

class EnaSubGui(Frame):
    # The AlleleGui class is an extension of Tkinter Frame. GUI elements go here.

    def selectAll(self, event):
        # I shouldn't need to write a select-All method but TK is kind of annoying.
        event.widget.tag_add("sel", "1.0", "end")

    def customPaste(self, event):
        # Writing a Paste method because TK doesn't operate intuitively. It doesn't replace the highlighted text.
        try:
            event.widget.delete("sel.first", "sel.last")
        except:
            pass
        event.widget.insert("insert", event.widget.clipboard_get())
        return "break"

    def __init__(self, root):
        # Initialize the GUI
        Frame.__init__(self, root)
        root.title("Create and Submit an EMBL-ENA Sequence Submission")
        self.parent = root
        
        # Assign the icon of this sub-window.
        assignIcon(self.parent)

        # Some basic functionality just doesn't work in TK. For example Ctrl-A, and pasting text.
        # I shouldn't need to encode this but here we are.
        root.bind_class("Text","<Control-a>", self.selectAll)
        root.bind_class("Text", "<<Paste>>", self.customPaste)
        
        # To define the exit behavior.  Save the input sequence text.
        self.parent.protocol('WM_DELETE_WINDOW', self.saveAndExit)

        button_opt = {'fill': BOTH, 'padx': 35, 'pady': 5}

        # Setup the Menu bar
        mainMenuBar = Menu(root)

        batchMenu = Menu(mainMenuBar, tearoff=0)
        batchMenu.add_command(label='Configure Batch Settings', command=self.chooseSubmissionOptions)
        batchMenu.add_command(label='Import Submission Data File', command=self.importSubmissions)
        batchMenu.add_command(label='Export Submission Data File', command=self.exportSubmissions)
        mainMenuBar.add_cascade(label='Batch Options', menu=batchMenu)

        helpMenu = Menu(mainMenuBar, tearoff=0)
        helpMenu.add_command(label='How To Use This tool', command=self.howToUse)
        helpMenu.add_command(label='Sample Submission Data', command=self.sampleSequence)
        mainMenuBar.add_cascade(label='Help', menu=helpMenu)

        root.config(menu=mainMenuBar)

        # A frame for the Instructions Label.
        self.instructionsFrame = Frame(self)  
        self.instructionText = StringVar()       
        self.instructionText.set('\nThis tool will generate a batch of HLA allele submissions for\n'
            + 'the EMBL-ENA nucleotide database.\n')
        Label(self.instructionsFrame, width=85, height=3, textvariable=self.instructionText).pack()
        self.instructionsFrame.pack(expand=False, fill='both')

        # Gather the submission batch to put the data on this frame.
        self.submissionIndex = 0 # This 0-based index indicates what submission within the batch we are currently looking at.
        self.submissionBatch = getConfigurationValue('submission_batch')
        logging.debug ('Just loaded the configuration. I have this many submissions in the batch:' + str(len(self.submissionBatch.submissionBatch)))
        # A frame to show what submission we are on in the batch.
        self.submIndexFrame = Frame(self)
        self.previousSubmissionButton = Button(self.submIndexFrame, text='<- Previous Submission', command=self.previousSubmission)
        self.previousSubmissionButton.grid(row=0, column=0)
        self.submissionIndexText = StringVar()
        # submission Index is zero-based. Add 1 for it to make sense to humans.
        self.submissionIndexText.set('.........')
        self.SubmissionIndexLabel = Label(self.submIndexFrame, width=25, height=3, textvariable=self.submissionIndexText)
        self.SubmissionIndexLabel.grid(row=0, column=1)
        self.nextSubmissionButton = Button(self.submIndexFrame, text='Next Submission ->', command=self.nextSubmission)
        self.nextSubmissionButton.grid(row=0, column=2)
        self.submIndexFrame.pack()

        self.addSubmissionFrame = Frame(self)
        self.deleteSubmissionButton = Button(self.addSubmissionFrame, text='Delete Submission', command=self.deleteCurrentSubmission)
        self.deleteSubmissionButton.grid(row=0, column=0)
        self.insertSubmissionButton = Button(self.addSubmissionFrame, text='New Submission', command=self.newSubmission)
        self.insertSubmissionButton.grid(row=0, column=1)
        self.addSubmissionFrame.pack()

        self.setSubmissionButtonState() # Buttons must exist before i call this.

        # Create a frame for basic allele sequence information.
        self.sequenceIdentifierFrame = Frame(self)
        #Standard Inputs widths for the form elements
        formInputWidth = 30
        labelInputWidth = 30

        self.sampleIDInstrText = StringVar()
        self.sampleIDInstrText.set('Sample ID:')
        self.sampleIDinstrLabel = Label(self.sequenceIdentifierFrame, width=labelInputWidth, height=1, textvariable=self.sampleIDInstrText).grid(row=0, column=0)
        self.inputSampleID = StringVar()
        self.inputSampleIDEntry = Entry(self.sequenceIdentifierFrame, width=formInputWidth, textvariable=self.inputSampleID)
        self.inputSampleIDEntry.grid(row=0, column=1)

        self.geneInstrStringVar = StringVar()
        self.geneInstrStringVar.set('Gene:')
        self.geneInstrLabel = Label(self.sequenceIdentifierFrame, width=labelInputWidth, height=1, textvariable=self.geneInstrStringVar).grid(row=1, column=0)
        self.inputGene = StringVar()
        self.inputGeneEntry = Entry(self.sequenceIdentifierFrame, width=formInputWidth, textvariable=self.inputGene)
        self.inputGeneEntry.grid(row=1, column=1)

        self.chooseClassIntVar = IntVar()
        self.chooseClassIntVar.set(1)
        self.classIRadio = Radiobutton(self.sequenceIdentifierFrame, text="HLA Class I ", variable=self.chooseClassIntVar, value=1)
        self.classIRadio.grid(row=2, column=0)
        self.classIIRadio = Radiobutton(self.sequenceIdentifierFrame, text="HLA Class II", variable=self.chooseClassIntVar, value=2)
        self.classIIRadio.grid(row=2, column=1)

        self.alleleInstrText = StringVar()
        self.alleleInstrText.set('Allele Local Name:')
        self.alleleInstrLabel = Label(self.sequenceIdentifierFrame, width=labelInputWidth, height=1, textvariable=self.alleleInstrText).grid(row=3, column=0)
        self.inputAllele = StringVar()
        self.inputAlleleEntry = Entry(self.sequenceIdentifierFrame, width=formInputWidth, textvariable=self.inputAllele)
        self.inputAlleleEntry.grid(row=3, column=1)

        self.sequenceIdentifierFrame.pack()

        # Create a frame for the input widget and scrollbars.
        self.featureInputFrame = Frame(self)
        self.featureInstrText = StringVar()
        self.featureInstrText.set('Annotated Sequence:')
        self.featureInstrLabel = Label(self.featureInputFrame, width=80, height=1, textvariable=self.featureInstrText).pack()
        self.featureInputXScrollbar = Scrollbar(self.featureInputFrame, orient=HORIZONTAL)
        self.featureInputXScrollbar.pack(side=BOTTOM, fill=X)
        self.featureInputYScrollbar = Scrollbar(self.featureInputFrame)
        self.featureInputYScrollbar.pack(side=RIGHT, fill=Y)
        self.featureInputGuiObject = Text(
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
        self.submButtonFrame = Frame(self)

        #self.submissionOptionsButton = Button(self.submButtonFrame, text='1) Submission Options', command=self.chooseSubmissionOptions)
        #self.submissionOptionsButton.grid(row=0, column=0)
        self.annotateFeaturesButton = Button(self.submButtonFrame, text='2) Annotate Exons & Introns' , command=self.annotateInputSequence)
        self.annotateFeaturesButton.grid(row=0, column=1)
        #self.generateSubmissionButton = Button(self.submButtonFrame, text='3) Generate an ENA-ENA submission', command=self.constructSubmission)
        #self.generateSubmissionButton.grid(row=0, column=2)
        self.submButtonFrame.pack()

        # Output interface is contained on a frame.
        self.submOutputFrame = Frame(self)
        self.outputENASubmission = StringVar()
        self.outputENASubmission.set('Allele Submission Preview:')
        self.outputENALabel = Label(self.submOutputFrame, width=80, height=1, textvariable=self.outputENASubmission).pack()
        self.submOutputXScrollbar = Scrollbar(self.submOutputFrame, orient=HORIZONTAL)
        self.submOutputXScrollbar.pack(side=BOTTOM, fill=X)
        self.submOutputYScrollbar = Scrollbar(self.submOutputFrame)
        self.submOutputYScrollbar.pack(side=RIGHT, fill=Y)
        self.submOutputGuiObject = Text(
            self.submOutputFrame, width=80, height=8, wrap=NONE
            , xscrollcommand=self.submOutputXScrollbar.set
            , yscrollcommand=self.submOutputYScrollbar.set
            , state=DISABLED
        )
        self.submOutputXScrollbar.config(command=self.submOutputGuiObject.xview)
        self.submOutputYScrollbar.config(command=self.submOutputGuiObject.yview)
        self.submOutputGuiObject.pack(expand=True, fill='both') 
        self.submOutputFrame.pack(expand=True, fill='both')


        self.uploadSubmissionFrame = Frame(self)
        # TODO: Think about: Should I export a zip file with all the submissions?
        # TODO: Or else a .csv file with the submision data?
        # TODO: Duh! All of it! Submission files and data in a .zip.
        # TODO: In that case, remove the export button here. I put it above, in the GUI.
        #self.saveSubmissionButton = Button(self.uploadSubmissionFrame, text='4) Save Submission to My Computer', command=self.saveSubmissionFile)
        #self.saveSubmissionButton.pack(**button_opt)
        self.uploadButton = Button(self.uploadSubmissionFrame, text='5) Upload Submission to EMBL-ENA', command=self.uploadSubmission)
        self.uploadButton.pack(**button_opt)
        self.exitButton = Button(self.uploadSubmissionFrame, text='Exit', command=self.saveAndExit)
        self.exitButton.pack(**button_opt)
        self.uploadSubmissionFrame.pack()
        
        self.pack(expand=True, fill='both')

        self.loadCurrentSubmission()

    def loadCurrentSubmission(self):
        currentSubmission = self.submissionBatch.submissionBatch[self.submissionIndex]

        newSubmissionIndexText = 'Submission ' + str(self.submissionIndex + 1) + ' / ' + str(len(self.submissionBatch.submissionBatch))
        self.submissionIndexText.set(newSubmissionIndexText)

        self.setSubmissionButtonState()

        self.inputSampleID.set('' if currentSubmission.cellId is None else currentSubmission.cellId)
        self.inputGene.set('' if currentSubmission.submittedAllele.geneLocus is None else currentSubmission.submittedAllele.geneLocus)
        self.inputAllele.set('' if currentSubmission.localAlleleName is None else currentSubmission.localAlleleName)
        self.chooseClassIntVar.set(int(currentSubmission.submittedAllele.hlaClass) if currentSubmission.submittedAllele.hlaClass is not None else 1) # I think this works? 1 or 2 is always stored. 1 or 2 are the potential values of the intvar.
        self.overwriteSequenceText(' ')
        self.overwriteSubmissionText(' ')
        self.overwriteSequenceText(currentSubmission.submittedAllele.getAnnotatedSequence(includeLineBreaks=True))
        self.constructSubmission()

    def saveCurrentSubmission(self):
        currentSubmission = self.submissionBatch.submissionBatch[self.submissionIndex]

        currentSubmission.submittedAllele.rawSequence = self.featureInputGuiObject.get('1.0', 'end')
        currentSubmission.submittedAllele.identifyFeaturesFromFormattedSequence()
        currentSubmission.localAlleleName = self.inputAllele.get()
        currentSubmission.submittedAllele.geneLocus = self.inputGene.get()
        currentSubmission.cellId = self.inputSampleID.get()
        currentSubmission.submittedAllele.hlaClass = str(self.chooseClassIntVar.get())# I think this works? 1 or 2 is always stored. 1 or 2 are the potential values of the intvar.

        self.submissionBatch.submissionBatch[self.submissionIndex] = currentSubmission

    def importSubmissions(self):
        # Popup warning message: This will ADD the imported submissions to the batch.
        # Figure out a way to either replace or add the submission files.

        logging.error('Have not implemented import submission files yet.')

    def exportSubmissions(self):
        logging.error('Have not implemented export submission files yet.')

    def deleteCurrentSubmission(self):
        logging.debug('deleteCurrentSubmission pressed')
        self.saveCurrentSubmission()

        if(len(self.submissionBatch.submissionBatch) == 1):
            showInfoBox('Cannot delete last submission.','You cannot delete the last remaining submission in the batch.')
        else:
            del self.submissionBatch.submissionBatch[self.submissionIndex]

            # If that was the rightmost in the batch, we need to reduce the index
            if ((self.submissionIndex) >= len(self.submissionBatch.submissionBatch)):
                self.submissionIndex = self.submissionIndex - 1

            self.loadCurrentSubmission()

    def newSubmission(self):
        logging.debug('newSubmission pressed')
        self.saveCurrentSubmission()
        # I'm adding a submission AFTER the current one. This feels intuitive to me.
        self.submissionBatch.submissionBatch.insert(self.submissionIndex + 1,AlleleSubmission())
        self.submissionIndex += 1
        self.loadCurrentSubmission()

    def previousSubmission(self):
        logging.debug('previousSubmission pressed')
        self.saveCurrentSubmission()
        if(self.submissionIndex <= 0):
            # Button should be disabled, but we don't want to walk left if we're at minimum.
            pass
        else:
            self.submissionIndex = self.submissionIndex - 1
            self.loadCurrentSubmission()

    def nextSubmission(self):
        logging.debug('nextSubmission pressed')
        self.saveCurrentSubmission()
        if((self.submissionIndex + 1) >= len(self.submissionBatch.submissionBatch)):
            # Button should be disabled, but we don't want to walk right if we're at maximum.
            pass
        else:
            self.submissionIndex = self.submissionIndex + 1
            self.loadCurrentSubmission()

    def chooseSubmissionOptions(self):
        logging.info ('Opening the EMBL-ENA Submission Options Dialog')
        
        self.disableGUI()
        
        enaOptionsRoot = Toplevel()
        enaOptionsRoot.bind("<Destroy>", self.enableGUI)
        EnaSubOptionsForm(enaOptionsRoot).pack()
        
        # Set the X and the Y Position of the options window, so it is nearby.  
        enaOptionsRoot.update()
        windowXpos = str(self.parent.winfo_geometry().split('+')[1])
        windowYpos = str(self.parent.winfo_geometry().split('+')[2])
        newGeometry = (str(enaOptionsRoot.winfo_width()) + 'x'
            + str(enaOptionsRoot.winfo_height()) + '+'
            + str(windowXpos) + '+' 
            + str(windowYpos))
        enaOptionsRoot.geometry(newGeometry)

        enaOptionsRoot.mainloop()
        
    def sampleSequence(self):
        logging.debug('sampleSequence pressed')
        logging.error('Have not implemented sample sequence yet.')

        # TODO: Old code. I need to create the batch submission here.
        # self.featureInputGuiObject.insert('1.0', 'aag\nCGTCGT\nccg\nGGCTGA\naat')
        #
        # # Clear the password, keep the username
        # #assignConfigurationValue('ena_username','')
        # assignConfigurationValue('ena_password','')
        #
        # assignConfigurationValue('sample_id', 'sample_12345')
        # assignConfigurationValue('gene','HLA-C')
        # assignConfigurationValue('class','1')
        # assignConfigurationValue("allele_name",'Allele*01:02MstNew.6')
        #
        # assignConfigurationValue('study_accession','PRJEB12345')
        #
        # assignConfigurationValue('choose_project','2')
        #
        # assignConfigurationValue('study_identifier','HLA_Analysis_Project')
        # assignConfigurationValue('study_short_title','HLA Typing for Cancer Research.')
        # assignConfigurationValue('study_abstract','An abstract is a more in-depth description of the nature of the research project.')
        #
        # assignConfigurationValue('analysis_alias','unique_HLA_analysis_alias')
        # assignConfigurationValue('analysis_title','Novel HLA sequence from patient with Leukemia')
        # assignConfigurationValue('analysis_description','This is an HLA-A sequence from a patient. It was discovered that he has Leukemia, so we decided to sequence his HLA.')
        #
        # self.constructSubmission()

    def howToUse(self):
        # This method should popup some instruction text in a wee window.
        # This should be explicit on how to use the tool.
        logging.error('howToUse() is probably outdated. Check if it needs updating.')

        showInfoBox('How to use this tool',
            'This software is to be used to create an\n'
            + 'ENA-formatted submission document,\n'
            + 'which specifies a (novel) HLA allele.\n\n'       
                       
            + 'This tool requires you to submit a\n'
            + 'full length HLA allele, including\n'
            + '5\' and 3\' UTRs.\n\n'
     
            + 'To create & submit an EMBL-ENA submission:\n\n'
            + '1.) Paste a full-length HLA sequence in\n'
            + 'the Annotated Sequence text area.\n'
            + '2.) Push [Submission Options] and provide\n'
            + 'the necessary sequence metadata.\n'
            + '3.) Push [Annotate Exons & Introns] to\n'
            + 'annotate your exons automatically.\n'
            + '4.) Push [Generate an EMBL-ENA submission]\n'
            + 'button to generate a submission.\n'
            + '5.) Push [Upload Submission to EMBL-ENA]\n'
            + 'to submit the sequence\n'
            + 'using ENA Webin REST interface\n\n'
            
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

    # TODO: Delete this method, use exportSubmissions instead.
    """
    def saveSubmissionFile(self):
        logging.error('Deprecated the save submission file method. But somehow I still called the method.')
        # Ask user for a output file location, and write the ENA submission to a file.
        # This takes the input from the output field, rather than generate a new submission.
        # So the user can edit the submission before or after saving it.

        
        self.dir_opt = options = {}       
        options['initialdir'] = expanduser("~")
        options['parent'] = self
        options['title'] = 'Specify your output file.'
        options['initialfile'] = 'ENA.HLA.Submission.txt'
        outputFileObject = filedialog.asksaveasfile(**self.dir_opt)
        submissionText = self.submOutputGuiObject.get('1.0', 'end')
        outputFileObject.write(submissionText)
        
        # TODO: Did I detect any exceptions? Maybe I don't have permission to write that file
        # I saw an error when i wrote to a network drive once. 
    """
     
    def annotateInputSequence(self):
        logging.debug('Annotating Input Sequence')
        try:
            self.disableGUI()
            self.update()   
            
            # Popup.  This uses NMDP BTM ACT tool to annotate sequences.
            if (messagebox.askyesno('Annotate Sequence?'
                , 'This will annotate your sequence using the\n'
                + 'NMDP: BeTheMatch Gene Feature\n'
                + 'Enumeration / Allele Calling Tool.\n\n'
                + 'Do you want to continue?'
                )):

                roughNucleotideSequence = collectAndValidateRoughSequence(self.featureInputGuiObject.get('1.0', 'end'))
                currentSubmission = self.submissionBatch.submissionBatch[self.submissionIndex]
                currentSubmission.submittedAllele.rawSequence = roughNucleotideSequence
                currentSubmission.submittedAllele.annotateSequenceUsingService(rawRequestURL=getConfigurationValue('nmdp_act_rest_address'))
                self.overwriteSequenceText(currentSubmission.submittedAllele.getAnnotatedSequence(includeLineBreaks=True))

            self.update()
            self.enableGUI()
            
        except Exception:
            showInfoBox('Error Annotating Input Sequence.'
                , str(exc_info()))
            self.update()
            self.enableGUI()
            raise

    def overwriteSequenceText(self, sequenceText=None):
        if(sequenceText is not None and len(sequenceText)>0):
            oldState = self.featureInputGuiObject.cget('state')
            self.featureInputGuiObject.config(state=NORMAL)
            self.featureInputGuiObject.delete('1.0', 'end')
            self.featureInputGuiObject.insert('1.0', sequenceText)
            self.featureInputGuiObject.config(state=oldState)
        else:
            logging.warning('Attempting to replace sequence text with an empty value.')

    def overwriteSubmissionText(self, submissionText=None):
        if (submissionText is not None and len(submissionText) > 0):
            oldState = self.submOutputGuiObject.cget('state')
            self.submOutputGuiObject.config(state=NORMAL)
            self.submOutputGuiObject.delete('1.0', 'end')
            self.submOutputGuiObject.insert('1.0', submissionText)
            self.submOutputGuiObject.config(state=oldState)
        else:
            logging.warning('Attempting to replace submission text with an empty value.')

            
    def uploadSubmission(self):
        logging.debug('Uploading Submission Batch')
        # TODO: Perform a Batch Submission
        #performFullSubmission(self.submOutputGuiObject.get('1.0', 'end') )
        logging.error('Upload Submission needs batch Submission functionality. Add it.')
        
    def constructSubmission(self):
        # Gather sequence information from the input elements, and generate a text ENA submission.
        logging.debug('Constructing Submission')
        try:
            currentSubmission = self.submissionBatch.submissionBatch[self.submissionIndex]

            # TODO: What happens when the sequence is not annotated yet? Probably need to remove this logic. Annoying.
            # if (isSequenceAlreadyAnnotated(roughNucleotideSequence)):
            #     annotatedSequence = roughNucleotideSequence
            #
            # else:
            #
            #     if (messagebox.askyesno('Auto - Annotate Exons?'
            #         , 'It looks like your sequence features have not been identified.\n' +
            #         'Would you like to annotate using NMDP: BeTheMatch\n' +
            #         'Gene Feature Enumeration Tool?')):
            #
            #         self.annotateInputSequence()
            #         annotatedSequence = collectRoughSequence(self.featureInputGuiObject)
            #     else:
            #         # You chose not to annotate.  Hope this works out for you.
            #         annotatedSequence = roughNucleotideSequence
                
            allGen = EnaSubGenerator()
            allGen.submission = currentSubmission
            allGen.submissionBatch = self.submissionBatch
            enaSubmissionText = allGen.buildENASubmission()
                        
            if (enaSubmissionText is None or len(enaSubmissionText) < 1):
                #showInfoBox('Empty submission text'
                #    ,'You are missing some required information.\n'
                #    + 'Try the \'Submission Options\' button.\n')
                logging.warning('Submission text is empty.')

                self.overwriteSubmissionText('')
            else:
                self.overwriteSubmissionText(enaSubmissionText)

        except KeyError:
            showInfoBox('Missing Submission Options'
                ,'You are missing some required information.\n'
                + 'Use the \'Submission Options\' button.\n'
                + 'Missing Data: ' + str(exc_info()))
           
        except HlaSequenceException:
            showInfoBox('I see a problem with Sequence Format.'
                , str(exc_info()))
           
        except Exception:
            showInfoBox('Error Constructing Submission.'
                , str(exc_info()))
            raise
            
    def saveAndExit(self):
        logging.debug('saveAndExit pressed')
        self.saveCurrentSubmission()
        assignConfigurationValue('submission_batch',self.submissionBatch)

        writeConfigurationFile()
        self.parent.destroy()
        
    def enableGUI(self, event=None):
        self.toggleGUI(True)  
        
    def disableGUI(self):
        self.toggleGUI(False)   
        
    def toggleGUI(self, isEnabled): 
        logging.debug('Toggling GUI Widgets:' + str(isEnabled))
         
        newState = (NORMAL if (isEnabled) else DISABLED)
        
        # Choosing the widgets individually, this makes the most sense I think.
        #self.submissionOptionsButton.config(state=newState)
        #self.generateSubmissionButton.config(state=newState)
        self.annotateFeaturesButton.config(state=newState)
        self.uploadButton.config(state=newState)
        self.exitButton.config(state=newState)
        self.previousSubmissionButton.config(state=newState)
        self.nextSubmissionButton.config(state=newState)
        self.deleteSubmissionButton.config(state=newState)
        self.insertSubmissionButton.config(state=newState)
        self.classIRadio.config(state=newState)
        self.classIIRadio.config(state=newState)
        self.inputSampleIDEntry.config(state=newState)
        self.inputGeneEntry.config(state=newState)
        self.inputAlleleEntry.config(state=newState)
        self.featureInputGuiObject.config(state=newState)

        # set the Submission next/previous buttons, only if the GUI is enabled.
        if (isEnabled):
            self.setSubmissionButtonState()

    def setSubmissionButtonState(self):
        # Enable or disable the next/previous buttons.
        if(self.submissionIndex <= 0):
            self.previousSubmissionButton.config(state=DISABLED)
        else:
            self.previousSubmissionButton.config(state=NORMAL)

        if((self.submissionIndex + 1) >= len(self.submissionBatch.submissionBatch)):
            self.nextSubmissionButton.config(state=DISABLED)
        else:
            self.nextSubmissionButton.config(state=NORMAL)

        if(len(self.submissionBatch.submissionBatch) <= 1):
            self.deleteSubmissionButton.config(state=DISABLED)
        else:
            self.deleteSubmissionButton.config(state=NORMAL)
