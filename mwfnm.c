#include "mwfnm.h"

/*
gcc -shared -fPIC libmwfnm.cpp -o mwfnm.so

loadlibrary('mwfnm.so')
calllib('mwfnm','testMe',2)

nData = 4;
arr = single(ones(nData,1));
pData = libpointer('singlePtr', arr);
calllib('mwfnm', 'useArray', pData, nData)

 */


int testMe(int k)
{
  return 2*k;
}

float useArray(const float* pData, const int nData)
{
  float retval = 0.0f;
  int i;
  for (i = 0 ; i < nData ; i++) {
    retval += pData[i];
  }
  retval /= (float) nData;
  return retval;
}
