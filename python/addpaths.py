import sys
import os
import glob

# Improved version
filedir = os.path.dirname(os.path.realpath(__file__))
cwd     = os.getcwd()

# Go one levels up (if possible)
os.chdir(filedir)
iDir = 1
while (iDir > 0):
    os.chdir('..')
    iDir = iDir - 1;

matches = []

# TODO: Python 3.0 use '**/fnm/Release', recursive=True):
for i in range(3):
  if os.name == 'nt':
    paths = [i*'*/' + 'fnm',i*'*/' + 'fnm', i*'*/' + 'fnm']
  else:
    paths = [i*'*/' + 'fnm', i*'*/' + 'release/fnm', i*'*/' + 'debug/fnm']
  paths.append(i*'*/' + 'python')
  for path in paths:
    entries = glob.glob(path)
    if len(entries) > 0:
      path = os.path.join(os.getcwd(),entries[0])
      matches.append(path)

# Field (build in source)
matches.append(os.path.join(os.getcwd(),'f2'))

# Add installed version (last resort)
paths =['lib','python','bin']
for path in paths:
  entries = glob.glob(path)
  if len(entries) > 0:
    path = os.path.join(os.getcwd(),entries[0])
    matches.append(path)

for path in matches:
  sys.path.append( path )

# HACK
sys.path.append('/home/jmh/git/bft4/build/bftx')

os.chdir(cwd)
