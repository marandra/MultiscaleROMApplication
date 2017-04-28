#!/usr/bin/python
'''
Gappy selection of integration points for the Second Reduction HP-ROM
'''

from lsqnonneg import *
import scipy.io as sio
import time

factorLEQ=1.0
tol = 1e-10

nGP = 40


k = numpy.array([[0.1210, 0.2319, 0.4398, 0.9342, 0.1370],
                     [0.4508, 0.2393, 0.3400, 0.2644, 0.8188],
                     [0.7159, 0.0498, 0.3142, 0.1603, 0.4302],
                     [0.8928, 0.0784, 0.3651, 0.8729, 0.8903],
                     [0.2731, 0.6408, 0.3932, 0.2379, 0.7349],
                     [0.2548, 0.1909, 0.5915, 0.6458, 0.6873],
                     [0.8656, 0.8439, 0.1197, 0.9669, 0.3461],
                     [0.2324, 0.1739, 0.0381, 0.6649, 0.1660],
                     [0.8049, 0.1708, 0.4586, 0.8704, 0.1556],
                     [0.9084, 0.9943, 0.8699, 0.0099, 0.1911]])

GaussWeights = numpy.array([0.4225, 0.8560, 0.4902, 0.8159, 0.4608, 0.4574, 0.4507, 0.4122, 0.9016, 0.0056])

'''
k = numpy.random.rand(400000,500)
GaussWeights = numpy.random.rand(400000)
nGP = 300
'''

#mat = scipy.io.loadmat('file.mat')


# ==========================================================================================

# Compute Matrix J and RHS vector
[J,b,INTexact] = ComputeJandb(k,GaussWeights,factorLEQ)

#J=sio.loadmat('J.mat')['J']
#print J

#W=sio.loadmat('W.mat')['W']
#print W

M=len(GaussWeights)
y = numpy.arange(M)

# Resudual vector, initial guess
r = b

# Number of iterations
it = 0
mPOS=0

z = []

Jnorm = numpy.sqrt(sum(numpy.multiply(J,J),0))

print "Point Selection Algorithm"

while numpy.linalg.norm(r)/numpy.linalg.norm(b)>tol and mPOS<=nGP:

   # 1. Compute new point
   ObjFun = numpy.dot((J[:,y]).T, r)
   div = numpy.multiply(Jnorm[y],numpy.linalg.norm(r))

   ObjFun = numpy.divide(ObjFun,div)

   s = ObjFun.argmax()
   t =y[s]

   # 2. Move i from set y to set z
   z = (numpy.append(z,t)).astype(int)
   y=numpy.delete(y,s)

   #   #solving LS conventional problem
   x = numpy.linalg.lstsq(J[:,z],b)[0]

   if any(x<0):
      print "Solving NNLS"

      # 3. Determime alpha for solving a NNLS
      [x, resnorm, residual] = lsqnonneg(J[:,z], b)

   # 3. Determime alpha for solving a NNLS
   #[x, resnorm, residual] = lsqnonneg(J[:,z], b)

   # 4. Update the residual 
   r = b-numpy.dot(J[:,z],x)

   # 5. Update mPOS and k
   mPOS = len(numpy.where(x>0)[0])

   # Iteration counter
   it = it + 1

   print "k =", it, "--- mPOS =",mPOS, "--- error (%) = ", numpy.linalg.norm(r)/numpy.linalg.norm(b)*100
   

# 6. Postprocess of points - neglecting null weights
INDzero=numpy.where(x==0)[0]

if any(INDzero):
   z=numpy.delete(z,INDzero)   

w = numpy.multiply(x,numpy.sqrt(GaussWeights[z]))

print "Reduced Weights"
print w

print "GP's index"
print z
