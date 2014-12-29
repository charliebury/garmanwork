# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 16:58:21 2014

@author: lina2532
"""

import smtplib
from datetime import datetime

def scriptstart():
    startTime = datetime.now()
    return startTime


def scriptdone(startTime):
    script_duration = datetime.now() - startTime
    fromaddr = 'charliebury3@gmail.com'
    toaddrs  = 'csbury@me.com'

    # Send email
    senddate=datetime.strftime(datetime.now(), '%Y-%m-%d')
    subject="Your job has completed"
    m="Date: %s\r\nFrom: %s\r\nTo: %s\r\nSubject: %s\r\nX-Mailer: My-Mail\r\n\r\n" % (senddate, fromaddr, toaddrs, subject)
    msg='''
    Job runtime: '''+str(script_duration)
    
    # Credentials (if needed)
    username = 'charliebury3@gmail.com'
    password = 'Maujei123'

    # The actual mail send
    server = smtplib.SMTP('smtp.gmail.com:587')
    server.starttls()
    server.login(username,password)
    server.sendmail(fromaddr, toaddrs, m+msg)
    server.quit()