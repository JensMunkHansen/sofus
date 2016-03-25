void tripleIndex(int singleIndex, int indexArray[3], int
		 sizeI, int sizeJ, bool baseOfReturnIndex){
  /*The function trippleindex calculates three zero based
    (baseOfReturnIndex=0) or 1-based (baseOfReturnIndex
    =1) indices from a C++ MEX 0-based singleIndex, and
    the 1st and 2nd array lengths.*/
  //Foundational eq. is: singleIndex= i + sizeI*j + sizeI*sizeJ*k;
  int IxJ = sizeI * sizeJ;
  int rem = singleIndex % IxJ; //i.e. rem=sizeI*j + i;
  indexArray[0]=( rem % sizeI ); //This is the final i value (1st index)
  indexArray[1]=( ( rem - indexArray[0] ) / sizeI ); // This is the final j value (2nd index)
  indexArray[2]=( ( singleIndex - indexArray[0] -
		    indexArray[1] * sizeI ) / IxJ ); //This is the final k value (3rd index)
  if (baseOfReturnIndex){
    indexArray[0]++; indexArray[1]++;indexArray[2]++;
  }; //Modification of output if 1-based return indexing (MATLAB index) is specified (defualt return value is 0-based (C++ index))
  return;
};
