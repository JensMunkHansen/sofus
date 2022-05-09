import os

def nargout(*args):
  """
  nargout function similar to what Matlab has

  Usage:
    def fun()
      # do something
      # ...
      return nargout(3)

    r = fun() # 3 is assigned to r
    fun()     # nothing is assigned 
  """
  import traceback
  callInfo = traceback.extract_stack()
  if os.name == 'posix':
    callLine = str(callInfo[-3].line)
    split_equal = callLine.split('=')
    split_comma = split_equal[0].split(',')
    num = len(split_comma)
  else:
    tmp = traceback.format_list(callInfo)
    callLine = tmp[-3].split('\n')[-2]
    split_equal = callLine.split('=')
    if type(split_equal) == list:
      if len(split_equal) > 1:
        split_comma = split_equal[0].split(',')
        num = len(split_comma)
        if num > 0:
          return args[0:num] if num > 1 else args[0]

if __name__ == '__main__':
  import numpy as np
  def fun(nparray):
    return nargout(np.mean(nparray),
                   np.std(nparray))

  arr = np.array([3, 4, 5])
  mean = fun(arr)
  mean, std = fun(arr)
  print(mean)
  print(std)
