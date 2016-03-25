import numpy as np
import addpaths
import swig_fnm as lib
import sys

n = 128

img = np.zeros((n,n),dtype=np.uint8)
#img[68:70,63:65] = 1
img[2:4,2:4] = 1

uimg = lib.uint8ArrayClass(n*n)

img = img.flatten()

for i in range(n**2):
 lib.uint8Array_setitem(uimg,i,int(img[i]))

block = np.zeros((4,4),dtype=np.uint8)
block[1:3,1:3] = 1
blockf = block.flatten()

ublock = lib.uint8ArrayClass(4*4*4)
for j in range(4):
  for i in range(4*4):
    lib.uint8Array_setitem(ublock,i+j*4*4,int(blockf[i]))

block4 = np.zeros((4,4,4),dtype=np.uint8)

block4[:,0,:] = np.roll(block,-1,1)
block4[:,1,:] = np.roll(block,-1,0)
block4[:,2,:] = np.roll(block,1,0)
block4[:,3,:] = np.roll(block,1,1)

block4[:,0,:] 
for i in range(4): # row
  for j in range(4): # Block index
    for k in range(4): # col
      lib.uint8Array_setitem(ublock,k + j*4 + i*4*4,int(block4[i,j,k]))

      
match = lib.int32ArrayClass(8)
 
lib.blockMatch4x4(uimg,128,ublock,4,match,128,128)

print('match: %d, %d' % (match[0],match[1]))

# Figure out layout of block
lib.blockMatch4x4SSE2(uimg,128,ublock,16,match,128,128)

sys.stdout.write('match: ')
for i in range(8):
    sys.stdout.write('%d ' % match[i])
sys.stdout.write('\n')

lib.blockMatch4x4SSE4(uimg,128,ublock,16,match,128,128)

sys.stdout.write('match: ')
for i in range(8):
    sys.stdout.write('%d ' % match[i])
sys.stdout.write('\n')
