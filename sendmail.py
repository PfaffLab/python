#!/usr/bin/env python
#==============================================================================
# sendmail.py
#
# Shawn Driscoll
# 20131125
#
# Gene Expression Laboratory, Pfaff
# Salk Institute for Biological Studies
#
# Snagged from an ubuntu forums post about using gmail from bash and
# sending file attachments. i edited it so it will send mail using Salk's 
# SMTP so I can have this send me emails when analysis completes or whatever.
#
# link: http://ubuntuforums.org/showthread.php?t=1472520
#
#==============================================================================

# useage: python send.py "receiver@gmail.com" "My Topic" "My text" "My attachment.jpg"

import os, re
import sys
import smtplib
import argparse
 
#from email.mime.image import MIMEImage
from email import encoders
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.MIMEText import MIMEText

#==============================================================================
# main
#==============================================================================

def main(args):

	# variables
    attach_file = False

	# check input file
    if args.attachment is not None:
        attach_file = True
    	if not file_exists(args.attachment):
    		sys.stderr.write("[main] Error: input file doesn't exist (%s)\n" % args.infile)
    		return 1

    recipient = args.to
    subject = args.subject
    message = args.message


    msg = MIMEMultipart()
    msg['Subject'] = subject
    msg['To'] = recipient
    msg['From'] = sender
    
    
    part = MIMEText('text', "plain")
    part.set_payload(message)
    msg.attach(part)
    
    session = smtplib.SMTP(SMTP_SERVER, SMTP_PORT)
 
    session.ehlo()
#    session.starttls()
    session.ehlo
    
    session.login(sender, password)

    if attach_file:
        fp = open(args.attachment, 'rb')
        msgq = MIMEBase('audio', 'audio')
        msgq.set_payload(fp.read())
        fp.close()
        # Encode the payload using Base64
        encoders.encode_base64(msgq)
        # Set the filename parameter
        filename=args.attachment
        msgq.add_header('Content-Disposition', 'attachment', filename=filename)
        msg.attach(msgq)

    # Now send or store the message
    qwertyuiop = msg.as_string()

    session.sendmail(sender, recipient, qwertyuiop)
    
    session.quit()
    os.system('notify-send "Email sent"')

	#return 0


def file_exists(fname):
	try:
		fin = open(fname)
	except IOError as e:
		return False

	fin.close()
	return True

#==============================================================================
# entry point
#==============================================================================

parser = argparse.ArgumentParser(description="About.")
parser.add_argument('to', type=str, help="recipient")
parser.add_argument('subject', type=str, help="subject")
parser.add_argument('message', type=str, help="message")
parser.add_argument('-a', dest="attachment", type=str, default=None, help="file to attach (default: none)")

args = parser.parse_args()

#==-=-===----=-==-=-==-=-=-===---=-=-=-===--====---=--==--=---=--=---=--=-=---=
# globals
#==-=-===----=-==-=-==-=-=-===---=-=-=-===--====---=--==--=---=--=---=--=-=---=
 
#SMTP_SERVER = 'smtp.gmail.com'
#SMTP_PORT = 587

SMTP_SERVER = 'helix.salk.edu'
SMTP_PORT = 25 # or 825

sender = 'sdriscoll@salk.edu'
password = "S@!k12044"
#recipient = args.to
#subject = args.subject
#message = args.message

if __name__ == "__main__":
	sys.exit(main(args))
