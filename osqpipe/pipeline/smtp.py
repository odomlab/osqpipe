#
# Copyright 2018 Odom Lab, CRUK-CI, University of Cambridge
#
# This file is part of the osqpipe python package.
#
# The osqpipe python package is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# The osqpipe python package is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the osqpipe python package.  If not, see
# <http://www.gnu.org/licenses/>.

'''SMTP mail functions used for user and admin notifications.'''

from email.mime.text import MIMEText
import smtplib
from osqutil.config import Config
from ..models import User

CONFIG = Config()

def email_admins(subject, body):
  '''
  Send an email to the admin users as registered in our repository database.
  '''
  recips = set([ u.email for u in User.objects.filter(is_superuser=True) ])
  send_email(subject, body, recips)

def send_email(subject, body, recips, include_admins=False):
  '''
  Send an email to a list of recipients, optionally including the
  registered admin users.
  '''
  if include_admins:
    recips = set(recips)
    recips = recips.union([ u.email for u
                            in User.objects.filter(is_superuser=True) ])

  mail = MIMEText(body.encode('ascii', 'replace'))
  mail['Subject'] = subject.encode('ascii', 'replace')
  mail['From'] = CONFIG.smtp_sender
  mail['To'] = ",".join(recips)

  conn = smtplib.SMTP(host=CONFIG.smtp_server)
  conn.sendmail(CONFIG.smtp_sender, recips, str(mail))
  conn.close()
