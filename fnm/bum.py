
nmax = 5
for n in range(1,nmax+1):
  print('n: %d' % n)
  for k in range(1,n):
    print('scale, uv: %d, i: %d, j: %d:' % (k,k,n-k))
  print('')
  for k in range(n+1,nmax+1):
    print('scale uv*: %d, i: %d, j: %d:' % (n,k,k-n))


#     1
#    1 1
#   1 2 1
#  1 3 3 1
# 1 4 6 4 1
# Local variables: #
# python-indent-offset: 2 #
# End: #

