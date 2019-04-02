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



from __future__ import print_function
import pickle
import os.path
import base64
import time
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from googleapiclient.errors import HttpError
from apiclient.http import MediaFileUpload
from email.mime.text import MIMEText

from saddlebags.AlleleSubCommon import getSaddlebagsDirectory

from os.path import expanduser, join

# If modifying these scopes, delete the file token.pickle.
# I need to send emails. I don't think I need to read emails, i prefer not to if i can help it.
# /auth/drive is full read write delete access. Not very comfortable with that but I
# must be able to create files.
SCOPES = [
    'https://www.googleapis.com/auth/drive',
    'https://www.googleapis.com/auth/gmail.send'
]

# Commenting this, I'm not runnign this as a main method.
#if __name__ == '__main__':
#    main()



# Important setup step:
# https://developers.google.com/drive/api/v3/quickstart/python
# There is a link there to authorize and create a credential file.
# google around, there may be simpler ways to do this.
# https://developers.google.com/drive/api/v3/quickstart/python
# The links on that page give access to a "quickstart" project.
# That works...but it's dumb. I want to automatically create a Saddlebags project.
# I haven't found a way to directly link and create a Saddlebags API project.
# Can be done through the browser like this (a magic link:):
# https://console.developers.google.com/flows/enableapi?apiid=drive
# Ah okay. I can use the browser to activate the Project in Google.
# I think this is telling me to have a gmail account for the project owner (saddlebags.admin@gmail.com)
# Then I can distribute the client_secrets.json with the application.
# Then I can allow users to login and authorize saddlebags to access their stuff. (i think)
# https://cloud.google.com/docs/authentication/end-user
# Another tutorial: https://o7planning.org/en/11917/create-credentials-for-google-drive-api
# https://codelabs.developers.google.com/codelabs/gsuite-apis-intro/#7


def uploadZipToIpdHla(zipFileName):

    """Shows basic usage of the Drive v3 API.
    Prints the names and ids of the first 10 files the user has access to.
    """


    print ('Uploading a zip file to IMGT/HLA:' + str(zipFileName))


    homeTempDirectory = getSaddlebagsDirectory()
    zipFileFullPath = join(homeTempDirectory, zipFileName)
    #jsonFileLocation = join(homeTempDirectory, 'credentials.json')
    #tokenPickleLocation = join(homeTempDirectory, 'token.pickle')

    service = getGoogleDriveService()

    # Call the Drive v3 API
    # Get a list of 10 files.
    results = service.files().list(
        pageSize=10, fields="nextPageToken, files(id, name)").execute()
    items = results.get('files', [])

    # list the first 10 files.
    if not items:
        print('No files found.')
    else:
        print('Files:')
        for item in items:
            print(u'{0} ({1})'.format(item['name'], item['id']))


    # TODO: Check if the file with this name already exists.
    # Conversely, name the zip file with a timestamp. This is risky because
    # we might accidentally send dupe zip files.

    # TODO: Detect the "Saddlebags" folder name from the list of files I have access to.
    # Conversely - Ask IPD to setup the upload folder.
    # They will send the upload folder directory.
    # Dominic suggests, I set up the initial configuration for the IPD upload.
    # IPD will set up my account for bulk submissions.
    # IPD sets up a shared google drive folder, which they share with a google email that I provide.


    # TODO: Send an Email to ipd submissions that a file has been uploaded.
    # Use google interface to do this again.
    # put (Test) in the email subject.

    # This folder is shared with mumc.tissuetyping.
    #folder_id = '1haV3w8DYawFKo0V4-pwwXoai4JUFaqeT'

    # This folder is in the personal drive of mumc.transplantationimmunology
    folder_id = '1Yw6Skb_ZsQn2kCbYUjo_zTnjSf-jstxX'

    print('zip file full path:' + zipFileFullPath)
    # Lets try to upload the file
    file_metadata = {'name': zipFileName,
        'parents': [folder_id]}
    media = MediaFileUpload(zipFileFullPath,
                            mimetype='application/zip')
    file = service.files().create(body=file_metadata,
                                        media_body=media,
                                        fields='id').execute()
    print('File ID: %s' % file.get('id'))

    # TODO: Somehow export the gui stuff to the GUI. Popups shouldn't happen in here.
    # Popup a query. Do you want to send an email to ipd submissions?

    # If yes, send an email.
        # Subject = Saddlebags HLA Allele Submission
        # Subject = (TEST) Saddlebags HLA Allele Submission
        # Body = "IPD User (username) has uploade the file "HLASubmission.zip"
        # to the shared folder "Saddlebags" (with the folder id (X))

        # Body += "Since this is a test submission, this file is safe to remove."

    emailSubject = '(TEST) Saddlebags HLA Allele Submission'
    emailBody  = ('IPD user (username) has just uploaded the file ' + zipFileName + '\n' +
        'It was placed into the folder ' + 'Saddlebags' + ' with a google folder resource id of ' + folder_id)

    if(True):
        emailBody += '\nSince this is a test submission, you are free to remove the file.'


    # or should the sender be "me"...maybe it's silly to pass this information in.
    # I reckon I could hard-code the target email address. not sure what it is right now, ipdsubmissions@anthonynolan.org? I made that up.
    sendEmail(emailSubject, emailBody, 'ben.matern@mumc.nl', 'me')


def getService(serviceName, apiVersion):
    # For some reason, the drive API suggests 'v3' for the api version, but gmail wants 'v1'. Ok.


    creds = None
    # The file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    homeTempDirectory = getSaddlebagsDirectory()
    #jsonFileLocation = join(homeTempDirectory, 'credentials.json')
    jsonFileLocation = join(homeTempDirectory, 'client_secrets.json')
    tokenPickleLocation = join(homeTempDirectory, 'token.pickle')


    if os.path.exists(tokenPickleLocation):
        with open(tokenPickleLocation, 'rb') as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                jsonFileLocation, SCOPES)
            creds = flow.run_local_server()
        # Save the credentials for the next run
        with open(tokenPickleLocation, 'wb') as token:
            pickle.dump(creds, token)
    service = build(serviceName, apiVersion, credentials=creds)
    return service

def getGmailService():
    return getService('gmail', 'v1')

def getGoogleDriveService():
    return getService('drive', 'v3')

def sendEmail(subject, body, targetEmailAddress, senderEmailAddress):
#def create_message(sender, to, subject, message_text):
    """Create a message for an email.

    Args:
      sender: Email address of the sender.
      to: Email address of the receiver.
      subject: The subject of the email message.
      message_text: The text of the email message.

    Returns:
      An object containing a base64url encoded email object.
    """

    message = MIMEText(body)
    message['to'] = targetEmailAddress
    message['from'] = senderEmailAddress
    message['subject'] = subject

    # Some bugs with python 3 and encoding the message properly.
    # https://stackoverflow.com/questions/38633781/python-gmail-api-not-json-serializable
    #encodedMessage = {'raw': base64.urlsafe_b64encode(message.as_string().encode("utf-8"))}
    rawMessage = base64.urlsafe_b64encode(message.as_bytes())
    rawMessage = rawMessage.decode()

    encodedMessage = {'raw': rawMessage}



#def send_message(service, user_id, message):
    """Send an email message.
    
    Args:
    service: Authorized Gmail API service instance.
    user_id: User's email address. The special value "me"
    can be used to indicate the authenticated user.
    (i don't understand how to use "me")
    message: Message to be sent.
    
    Returns:
    Sent Message.
    """

    service = getGmailService()

    try:
        message = (service.users().messages().send(userId='me', body=encodedMessage).execute())
        print('Message Id: %s' % message['id'])
        return message
    except HttpError as err:
        print('An error occurred: %s' % err)
        if err.resp.status in [403, 500, 503]:
            sleep(5)
        else:
            raise