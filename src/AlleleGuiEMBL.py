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


import os
from os import makedirs
from os.path import expanduser, join, isdir

import datetime
import hashlib
import ftplib
import gzip
import shutil
#import pycurl
import requests

import Tkinter, Tkconstants, tkFileDialog, tkMessageBox
from Tkinter import *

from SubmissionGeneratorEMBL import SubmissionGeneratorEMBL
from AlleleGuiEMBLInputForm import AlleleGuiEMBLInputForm
from AlleleSubCommon import *
from AlleleSubmissionEMBLXml import *
#from HLAGene import HLAGene

# The AlleleGui class is an extension of Tkinter.  The GUI elements and interactions are specified in this class.
class AlleleGuiEMBL(Tkinter.Frame):

    # I shouldn't need to write a select-All method but TK is kind of annoying.
    def selectall(self, event):
        event.widget.tag_add("sel","1.0","end")
        
    # Initialize the GUI
    def __init__(self, root):
        Tkinter.Frame.__init__(self, root)
        root.title("Create and Submit an EMBL Sequence Submission")
        self.parent = root

        # Ctrl-A doesn't work by default in TK.  I guess I need to do it myself.
        root.bind_class("Text","<Control-a>", self.selectall)
        
        # To define the exit behavior.  Save the input sequence text.
        self.parent.protocol('WM_DELETE_WINDOW', self.saveAndExit)

        button_opt = {'fill': Tkconstants.BOTH, 'padx': 35, 'pady': 5}
        
        
        # A frame for the Instructions Label.
        self.instructionsFrame = Tkinter.Frame(self)  
        self.instructionText = Tkinter.StringVar()       
        self.instructionText.set('\nThis tool will generate an HLA allele submission for\n'
            + 'the EMBL / ENA nucleotide database.\n'
            + 'If you provide login credentials, you may automatically submit the sequence.\n'
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
        self.submissionOptionsButton = Tkinter.Button(self.submButtonFrame, text='Submission Options', command=self.chooseSubmissionOptions)
        self.submissionOptionsButton.grid(row=0, column=0)
        self.generateSubmissionButton = Tkinter.Button(self.submButtonFrame, text=unichr(8681) + ' Generate an EMBL submission ' + unichr(8681), command=self.constructSubmission)
        self.generateSubmissionButton.grid(row=0, column=1)
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
        self.uploadButton = Tkinter.Button(self.uploadSubmissionFrame, text='Upload Submission to EMBL', command=self.uploadSubmission)
        self.uploadButton.pack(**button_opt)
        self.saveSubmissionButton = Tkinter.Button(self.uploadSubmissionFrame, text='Save Submission to My Computer', command=self.saveSubmissionFile)
        self.saveSubmissionButton.pack(**button_opt)
        self.exitButton = Tkinter.Button(self.uploadSubmissionFrame, text='Exit', command=self.saveAndExit)
        self.exitButton.pack(**button_opt)
        self.uploadSubmissionFrame.pack()
        
        self.pack(expand=True, fill='both')
         
    def chooseSubmissionOptions(self):
        print ('Opening the EMBL Submission Options Dialog')
        
        self.disableGUI()
        
        emblOptionsRoot = Tkinter.Toplevel()
        emblOptionsRoot.bind("<Destroy>", self.enableGUI)
        AlleleGuiEMBLInputForm(emblOptionsRoot).pack()
        
        # Set the X and the Y Position of the options window, so it is nearby.  
        emblOptionsRoot.update()        
        windowXpos = str(self.parent.winfo_geometry().split('+')[1])
        windowYpos = str(self.parent.winfo_geometry().split('+')[2])
        newGeometry = (str(emblOptionsRoot.winfo_width()) + 'x' 
            + str(emblOptionsRoot.winfo_height()) + '+' 
            + str(windowXpos) + '+' 
            + str(windowYpos))
        emblOptionsRoot.geometry(newGeometry)

        emblOptionsRoot.mainloop()
        
        
    def writeMd5(self, inputFileName, outputFileName):
        hash_md5 = hashlib.md5()
        with open(inputFileName, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        hashValue= hash_md5.hexdigest()
        
        outputFile = createOutputFile(outputFileName)
        # The Ubuntu md5sum program seems to write a single checksum and filename with 2 spaces between
        # I don't know why 2 spaces, but I'll roll with it.
        outputFile.write(str(hashValue) + '  ' + str(split(inputFileName)[1]))
        outputFile.close()

    def uploadSubmission(self):
        print('Uploading Submission to EMBL')

        emblUsername = getConfigurationValue('embl_username')
        emblPassword = getConfigurationValue('embl_password')
        if(emblUsername is None 
            or len(emblUsername) < 1
            or emblPassword is None 
            or len(emblPassword) < 1):
            tkMessageBox.showinfo('Missing Login Credentials', 
                'You must provide EMBL username and password.\n'
                'Please use the "Submission Options" button.')
            return
           

        useTestServers = (int(getConfigurationValue('test_submission')) == 1)
        # Are you sure?
        if useTestServers:
            result = tkMessageBox.askquestion("Submit to TEST / DEMO environment", "You are about to submit a sequence to the\n\nTEST / DEMO EMBL environment.\n\nAre You Sure?", icon='warning')
        else:
            result = tkMessageBox.askquestion("Submit to LIVE / PROD environment", "You are about to submit a sequence to the\n\nLIVE / PROD EMBL environment.\n\nAre You Sure?", icon='warning')

        if result == 'yes':
            pass
        else:
            return
        
        # TODO: Existing project? Maybe I should check if the study/project exists, before I get started
        
        
        
        

        # Determine a working directory. Folder underneath executable called temp.
        try:
            workingDirectory = join(expanduser("~"), 'temp_upload_directory')
            print('I can work in this directory:' + workingDirectory)
            
            if not isdir(workingDirectory):
                print('Making Directory:' + workingDirectory)
                makedirs(workingDirectory)
        except Exception:
            print 'Cannot Initialize Working Directory'
            print sys.exc_info()[1]
            tkMessageBox.showinfo('Working Directory Error', 
                'Sorry, I failed to create this working directory:\n'
                + str(workingDirectory)
                + '\n and I cannot continue.\nMaybe this is a '
                + 'permissions issue, are these folders read only?\n' 
                +  str(sys.exc_info()[1]))
            return
        
        # Give my submission a filename. SOmething with a datetime stamp
        try:
            # This includes a "seconds" measure, should be pretty unique.
            dateTimeNow = '{:%Y_%m_%d_%H_%M_%S}'.format(datetime.datetime.now())
            submissionFileName = join(workingDirectory, 'HLA_Submission_' + dateTimeNow + '.txt')
      
        except Exception:
            print 'Cannot Assign File Name'
            print sys.exc_info()[1]
            tkMessageBox.showinfo('File Name Error', 
                'Sorry, I failed to create this file:\n'
                + str(submissionFileName)
                + '\n and I cannot continue.\n' 
                +  str(sys.exc_info()[1]))
            return
        
        # Write submission to a file
        try:
            submissionText = self.submOutputGuiObject.get('1.0', 'end')           
            
            outputFileObject = open(submissionFileName, 'w') 
            outputFileObject.write(submissionText) 
            outputFileObject.close()        
        
        except Exception:
            print 'Cannot Write Submission Flatfile'
            print sys.exc_info()[1]
            tkMessageBox.showinfo('Cannot Write Submission Flatfile', 
                'Sorry, I failed to create the submission file:\n'
                + str(submissionText)
                + '\n and I cannot continue.\nMaybe this is a '
                + 'permissions issue, are these folders read only?\n' 
                +  str(sys.exc_info()[1]))
            return
        
        # gzip the submission file.  Make a gz file.
        try:
            zippedFileName = submissionFileName + '.gz'
            
            with open(submissionFileName, 'rb') as fileIn, gzip.open(zippedFileName, 'wb') as fileOut:
                shutil.copyfileobj(fileIn, fileOut)
        
        except Exception:
            print 'Cannot Compress Submission File'
            print sys.exc_info()[1]
            tkMessageBox.showinfo('Cannot Compress Submission File', 
                'Sorry, I failed to compress the submission file:\n'
                + str(zippedFileName)
                + '\n and I cannot continue.\n' 
                +  str(sys.exc_info()[1]))
            return
        
        # Calculate an MD5SUM
        try:
            md5FileName = zippedFileName + '.md5'
            self.writeMd5(zippedFileName,md5FileName)
            
        except Exception:
            print 'Cannot Calculate MD5'
            print sys.exc_info()[1]
            tkMessageBox.showinfo('Cannot Calculate an Md5 checksum', 
                'Sorry, I failed to calculate an md5 checksum\nand I cannot continue.\n' 
                +  str(sys.exc_info()[1]))
            return

        # Use FTP  to send the file to EMBL
        try:
            if useTestServers:
                ftpServerAddress = getConfigurationValue('embl_ftp_upload_site_test')      
            else:
                ftpServerAddress = getConfigurationValue('embl_ftp_upload_site_prod')   
            
            #print ('attempting to open ftp connection')
            ftp = ftplib.FTP(ftpServerAddress)
            ftp.login(getConfigurationValue('embl_username'), getConfigurationValue('embl_password'))
            ftp.storbinary('STOR ' + '/' + split(zippedFileName)[1], open(zippedFileName, 'rb'), 1024)
            ftp.storbinary('STOR ' + '/' + split(md5FileName)[1], open(md5FileName, 'rb'), 1024)
            ftp.close()
            # is that it?  Easy.

        except Exception:
            print 'Cannot Upload to FTP site'
            print sys.exc_info()[1]
            tkMessageBox.showinfo('Cannot Upload to FTP site', 
                'Sorry, I failed to upload your submission files to the EMBL FTP site\nand I cannot continue.\n' 
                +  str(sys.exc_info()[1]))
            return
        
        
        # TODO: I Need a center_name.  This is based on the username but it goes into these submissions.
        # Maybe I can get that from REST
        
        
        # Handle the new project
        # effectively, study = project 
        # existing study = 1, new study = 2
        newProject = (getConfigurationValue('choose_project') == '2')
        if newProject:
            
            # Generate Project and Project Submission XML Files
            try:
                projectFileName = join(workingDirectory, 'project.xml')
                projectText = createProjectXML(projectFileName
                    , getConfigurationValue('study_name')
                    , getConfigurationValue('study_description')
                    , getConfigurationValue('study_abstract'))
                
                projectSubFileName = join(workingDirectory, 'project_submission.xml')
                projectSubmissionText = createProjectSubmissionXML('proj_sub_' + dateTimeNow
                    , projectSubFileName)
                
                print('I made this project text:\n' + projectText)
                print('I made this project submission text:\n' + projectSubmissionText)
                
            except Exception:
                print 'Cannot Create Project Submission XML'
                print sys.exc_info()[1]
                tkMessageBox.showinfo('Cannot Create Project Submission XML', 
                    'Sorry, I failed to create a project XML file\nand I cannot continue.\n' 
                    +  str(sys.exc_info()[1]))
                return
                        
            # Use REST to submit this project
            try:
                #"https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%20Webin-NNN%20PASSWORD"
                requestURL = (getConfigurationValue('embl_rest_address_test')
                    + '?auth=ENA%20'
                    + getConfigurationValue('embl_username')
                    + '%20'
                    + getConfigurationValue('embl_password')          
                              
                               
                    )
                
                
                # So i keep getting 415 errors. Trying to figure out what's up with REST.
                
                print ('URL=\n'+str(requestURL))
                #auth=HTTPBasicAuth(getConfigurationValue('embl_username'), getConfigurationValue('embl_password'))
                #requestUser=(getConfigurationValue('embl_username'), getConfigurationValue('embl_password'))
                requestData = {'SUBMISSION':projectSubmissionText
                    , 'PROJECT':projectText
                   # ,'content-type':'text/xml'
                    }
                
                requestHeaders = {'content-type':'application/xml'
                   # 'HTTP Name':'Content-Type'
            
                    }
                
               
                
                postResponse = requests.post(
                    requestURL
                    , params=requestData
                    , headers=requestHeaders

#                    , auth=requestUser
                    #, auth=(getConfigurationValue('embl_username'), getConfigurationValue('embl_password'))
                    )
                
                print ('the post response object:\n' + str(postResponse))
                
                print ('response status:\n' + str(postResponse.status_code))
                
                print('response text:' + str(postResponse.text))
                
                print('response content:' + str(postResponse.content))
                
                print('is response okay?:' + str(postResponse.status_code == requests.codes.ok))
                
                
                # headers have the login credentials?
                
                #  -k, --insecure      Allow connections to SSL sites without certs (H)
                #  -F, --form CONTENT  Specify HTTP multipart POST data (H)
                # those login credentials are passed plain text in the URL? Why?
                #curl -k -F "SUBMISSION=@sub.xml" -F "PROJECT=@project.xml" "https://www-test.ebi.ac.uk/ena/submit/drop-box/submit/?auth=ENA%20Webin-NNN%20PASSWORD"
                #curl -k -F "SUBMISSION=@"$submissionFileName -F "PROJECT=@"$projectFileName $testEnaSite"?auth=ENA%20"$userName"%20"$password > curlProjResults.xml
                
                
                #data = ['SUBMISSION': '@projectSubFileName'
                #    ,
                #    'tx': str(request.GET.get('tx')),
                #    'at': paypal_pdt_test
                #    ]

                #post = urllib.urlencode(data)

                #c = pycurl.Curl()
                #if(useTestServers):
                #c.setopt(pycurl.URL, getConfigurationValue('embl_rest_address_test'))
                #c.setopt(pycurl.HTTPHEADER, ['X-Postmark-Server-Token: API_TOKEN_HERE','Accept: application/json'])
                #c.setopt(pycurl.POST, 1)
                #c.setopt(pycurl.POSTFIELDS, data)
                
                #print ('about to perform curl:' + str(c))
                #c.perform()
                
            except Exception:
                print 'Cannot Submit Project XML'
                print sys.exc_info()[1]
                tkMessageBox.showinfo('Cannot Submit Project XML', 
                    'Sorry, I failed to submit the project XML file\nand I cannot continue.\n' 
                    +  str(sys.exc_info()[1]))
                return   
                            
            # Create a new project REST
            # Open Response, determine if success
            
            # If errors:
                # report Errors and give up
            # else:
                # store accession#        
            
            
        else: #(existing project)
            pass

        
        # else (existing project)
            # Store the project accession #
        
        # Generate XML Files for new sequence
        # REST the new XML files over to them.
        # Open response determine if success?
        # Gather Important Accession Numbers
        # Store accession number in our config file
        # Popup message with Results

        tkMessageBox.showinfo('Success uploading submission.', 
            'Everything worked fine.\n' 
            +'Congratulations. TODO: This message should have more information in it.')
        
        
        
        
    def sampleSequence(self):
        self.featureInputGuiObject.delete('1.0','end')
        self.featureInputGuiObject.insert('1.0', 'aag\nCGTCGT\nccg\nGGCTGA\naat')
        
        # Clear the password, keep the username
        #assignConfigurationValue('embl_username','')
        assignConfigurationValue('embl_password','')
        
        assignConfigurationValue('sample_id', 'Donor_12345')
        assignConfigurationValue('gene','HLA-C')
        assignConfigurationValue('class','1')
        assignConfigurationValue("allele_name",'Allele:01:02')

        assignConfigurationValue('study_accession','PRJEB12345')
                                 
        assignConfigurationValue('choose_project','2')
        
        assignConfigurationValue('study_name','HLA_Analysis_Project')
        assignConfigurationValue('study_description','Our laboratory performs HLA typing for Research')
        assignConfigurationValue('study_abstract','An abstract is a more in-depth description of the nature of the research project.')
        
        self.constructSubmission()
        
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

    # Ask user for a output file location, and write the EMBL submission to a file.
    # This takes the input from the output field, rather than generate a new submission.
    # So the user can edit the submission before or after saving it.
    def saveSubmissionFile(self):

        self.dir_opt = options = {}       
        options['initialdir'] = expanduser("~")
        options['parent'] = self
        options['title'] = 'Specify your output file.'
        options['initialfile'] = 'EMBL.HLA.Submission.txt'
        outputFileObject = tkFileDialog.asksaveasfile(**self.dir_opt)
        submissionText = self.submOutputGuiObject.get('1.0', 'end')
        outputFileObject.write(submissionText)
        
        # TODO: Did I detect any exceptions? Maybe I don't have permission to write that file
        # I saw an error when i wrote to a network drive once. 

        
    # Gather sequence information from the input elements, and generate a text EMBL submission.
    def constructSubmission(self):
        try:

            allGen = SubmissionGeneratorEMBL()
            roughFeatureSequence = self.featureInputGuiObject.get('1.0', 'end')

            allGen.inputSampleID = getConfigurationValue('sample_id')
            allGen.inputGene = getConfigurationValue('gene')
            allGen.inputAllele = getConfigurationValue('allele_name')
            allGen.inputClass = getConfigurationValue('class')            
            
            allGen.processInputSequence(roughFeatureSequence)
            enaSubmission = allGen.buildENASubmission()
                        
            if (enaSubmission is None or len(enaSubmission) < 1):
                tkMessageBox.showerror('Empty submission text'
                    ,'You are missing some required information.\n'
                    + 'Try the \'Submission Options\' button.\n')
                
                self.submOutputGuiObject.delete('1.0','end')    
                self.submOutputGuiObject.insert('1.0', '') 
            else:
                self.submOutputGuiObject.delete('1.0','end')    
                self.submOutputGuiObject.insert('1.0', enaSubmission) 
            
        except KeyError, e:
            tkMessageBox.showerror('Missing Submission Options'
                ,'You are missing some required information.\n'
                + 'Use the \'Submission Options\' button.\n'
                + 'Missing Data: ' + str(e))
            
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
        self.featureInputGuiObject.config(state=newState)
        self.submissionOptionsButton.config(state=newState)
        self.generateSubmissionButton.config(state=newState)
        self.submOutputGuiObject.config(state=newState)
        self.uploadButton.config(state=newState)
        self.saveSubmissionButton.config(state=newState)
        self.exitButton.config(state=newState)

     
        
        
        
