import sys
import os
import fnmatch

mydir = os.path.dirname(os.path.realpath(__file__))


# Go two levels up (os.path.split saturates to ['/', ''])
iDir = 1
while (iDir > 0):
    mydir = os.path.split(mydir)[0]
    iDir = iDir - 1;

matches = []
for root, dirnames, filenames in os.walk(mydir):
  for dirname in fnmatch.filter(dirnames, 'fnm'):
      matches.append(os.path.join(root, dirname))
  for dirname in fnmatch.filter(dirnames, 'Release'):
      matches.append(os.path.join(root, dirname))
  for dirname in fnmatch.filter(dirnames, 'Debug'):
      matches.append(os.path.join(root, dirname))
  for dirname in fnmatch.filter(dirnames, 'python'):
      matches.append(os.path.join(root, dirname))

for path in matches:
    sys.path.append( path )

