#!/usr/bin/env python

"""
A Python implementation of NNLS algorithm
 
References:
[1]  Lawson, C.L. and R.J. Hanson, Solving Least-Squares Problems, Prentice-Hall, Chapter 23, p. 161, 1974.
 
Contributed by Klaus Schuch (schuch@igi.tugraz.at)
based on MATLAB's lsqnonneg function
 
"""
import numpy

eps = 2.22e-16    # from matlab 

def lsqnonneg(C, d, x0=None, tol=None, itmax_factor=10):
    '''Linear least squares with nonnegativity constraints.
 
    (x, resnorm, residual) = lsqnonneg(C,d) returns the vector x that minimizes norm(d-C*x)
    subject to x >= 0, C and d must be real
    '''
 
    #eps = 2.22e-16    # from matlab
    def norm1(x):
        return abs(x).sum().max()
 
    def msize(x, dim):
        s = x.shape
        if dim >= len(s):
            return 1
        else:
            return s[dim]
 
    if tol is None: tol = 10*eps*norm1(C)*(max(C.shape)+1)

    C = numpy.asarray(C)
 
    (m,n) = C.shape
    P = numpy.zeros(n)
    Z = numpy.arange(1, n+1)
 
    if x0 is None: x=P
    else:
        if any(x0 < 0): x=P
        else: x=x0
 
    ZZ=Z
 
    resid = d - numpy.dot(C, x)
    w = numpy.dot(C.T, resid)
    
    outeriter=0; it=0
    itmax=itmax_factor*n
    exitflag=1

    # outer loop to put variables into set to hold positive coefficients
    while numpy.any(Z) and numpy.any(w[ZZ-1] > tol):
        
        outeriter += 1
        
        t = w[ZZ-1].argmax()
        t = ZZ[t]
         
        P[t-1]=t
        Z[t-1]=0

        PP = numpy.where(P <> 0)[0]+1
        ZZ = numpy.where(Z <> 0)[0]+1
        
        CP = numpy.zeros(C.shape)
        CP[:, PP-1] = C[:, PP-1]
        CP[:, ZZ-1] = numpy.zeros((m, msize(ZZ, 1)))

        z=numpy.dot(numpy.linalg.pinv(CP), d)

        z[ZZ-1] = numpy.zeros((msize(ZZ,1), msize(ZZ,0)))

        # inner loop to remove elements from the positve set which no longer belong
        while numpy.any(z[PP-1] <= tol):
            it += 1
 
            if it > itmax:
                max_error = z[PP-1].max()
                raise Exception('Exiting: Iteration count (=%d) exceeded\n Try raising the \
                                 tolerance tol. (max_error=%d)' % (it, max_error))
 
            QQ = numpy.where((z <= tol) & (P <> 0))[0]

            alpha = min(x[QQ]/(x[QQ] - z[QQ]))

            x = x + alpha*(z-x)
 
            ij = numpy.where((abs(x) < tol) & (P <> 0))[0]+1
            Z[ij-1] = ij
            P[ij-1] = numpy.zeros(max(ij.shape))
            PP = numpy.where(P <> 0)[0]+1
            ZZ = numpy.where(Z <> 0)[0]+1
 
            CP[:, PP-1] = C[:, PP-1]
            CP[:, ZZ-1] = numpy.zeros((m, msize(ZZ, 1)))
 
            z=numpy.dot(numpy.linalg.pinv(CP), d)
            z[ZZ-1] = numpy.zeros((msize(ZZ,1), msize(ZZ,0)))
 
        x = z
        resid = d - numpy.dot(C, x)
        w = numpy.dot(C.T, resid)

        #print w 

    return (x, sum(resid * resid), resid)


def ComputeJandb(Modes, weights, factorLEQ=None):

    if factorLEQ is None: factorLEQ = 1.0

    # Exact integral - numerical integration
    #INTexact = numpy.dot(k.T,l.T)
    INTexact = numpy.dot(Modes.T,weights)

    # Total microscale volume
    vol=numpy.sum(weights)
    sqrtVol =numpy.sqrt(weights)

    #print vol, sqrtVol

    # Matrix of modified modes (with zero integral)
    Xf = numpy.zeros(Modes.shape)

    # Loops over the initial modes
    Xf = numpy.subtract(Modes,INTexact/vol)
    Xf = numpy.multiply(Xf.T,sqrtVol).T
    #print Xf1.T
    #print sqrtVol
    #for i in range(Modes.shape[1]):
    #     Xf[:,i] = Modes[:,i] - INTexact[i]/vol
    #     for j in range(Modes.shape[0]):
             #print Xf[:,i].shape
             #print j
             #print Xf[j,i]
    #         Xf[j,i] = Xf[j,i]*sqrtVol[j]
             #print Xf

    #print numpy.subtract(Xf1,Xf)
    #print Xf

    # Singular Value Decomposition
    [Lambda,SValues,VValues] = numpy.linalg.svd(Xf,full_matrices=False)

    # fixed tolerance to define the reduced modified set of modes
    tol = numpy.max(Modes.shape)*eps*numpy.max(SValues)
    RankXf = sum(i>tol for i in SValues)

    Lambda = Lambda[:,0:RankXf]

    J=Lambda.T
    Jw=factorLEQ*sqrtVol/numpy.sqrt(vol)

    # Adding last row related with the sqrt of gauss integration weigths
    J = numpy.vstack([J, Jw.T])

    # Initializing the RHS vector for the optimization problem
    b=numpy.append(numpy.zeros(Lambda.T.shape[0]),factorLEQ*numpy.sqrt(vol))

    return (J,b,INTexact)


 
# Unittest
if __name__=='__main__':
#    C = numpy.array([[0.0372, 0.2869],
#                     [0.6861, 0.7071],
#                     [0.6233, 0.6245],
#                     [0.6344, 0.6170]])
 
#    C1 = numpy.array([[0.0372, 0.2869, 0.4],
#                      [0.6861, 0.7071, 0.3],
#                      [0.6233, 0.6245, 0.1],
#                      [0.6344, 0.6170, 0.5]])
 
#    C2 = numpy.array([[0.0372, 0.2869, 0.4],
#                      [0.6861, 0.7071,-0.3],
#                      [0.6233, 0.6245,-0.1],
#                      [0.6344, 0.6170, 0.5]])
 
#    d = numpy.array([0.8587, 0.1781, 0.0747, 0.8405])
 
#    [x, resnorm, residual] = lsqnonneg(C, d)
#    dres = abs(resnorm - 0.8315)          # compare with matlab result
#    print 'ok, diff:', dres
#    if dres > 0.001:
#        raise Exeption('Error')
 
#    [x, resnorm, residual] = lsqnonneg(C1, d)
#    dres = abs(resnorm - 0.1477)          # compare with matlab result
#    print 'ok, diff:', dres
#    if dres > 0.01:
#        raise Exeption('Error')
 
#    [x, resnorm, residual] = lsqnonneg(C2, d)
#    dres = abs(resnorm - 0.1027)          # compare with matlab result
#    print 'ok, diff:', dres
#    if dres > 0.01:
#        raise Exeption('Error')
 
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


    k2= numpy.array([[0.814723686393179,         0.157613081677548,         0.655740699156587,         0.706046088019609,  0.438744359656398],
         [0.905791937075619,         0.970592781760616,        0.0357116785741896,        0.0318328463774207,         0.381558457093008],
         [0.126986816293506,         0.957166948242946,         0.849129305868777,          0.27692298496089,         0.765516788149002],
         [0.913375856139019,         0.485375648722841,         0.933993247757551,        0.0461713906311539,         0.795199901137063],
          [0.63235924622541,           0.8002804688888,         0.678735154857773,        0.0971317812358475,         0.186872604554379],
        [0.0975404049994095,         0.141886338627215,         0.757740130578333,         0.823457828327293,         0.489764395788231],
         [0.278498218867048,         0.421761282626275,         0.743132468124916,         0.694828622975817,         0.445586200710899],
         [0.546881519204984,         0.915735525189067,         0.392227019534168,         0.317099480060861,         0.646313010111265],
         [0.957506835434298,         0.792207329559554,         0.655477890177557,         0.950222048838355,         0.709364830858073],
         [0.964888535199277,         0.959492426392903,         0.171186687811562,        0.0344460805029088,         0.754686681982361]])



 
#    k1 = k-0.5
 
    l = numpy.array([0.4225, 0.8560, 0.4902, 0.8159, 0.4608, 0.4574, 0.4507, 0.4122, 0.9016, 0.0056])
 
    l2 = numpy.array([0.276025076998578, 0.679702676853675, 0.655098003973841, 0.162611735194631, 0.118997681558377, 0.498364051982143,  0.959743958516081,  0.340385726666133,  0.585267750979777,  0.223811939491137])    

    [x, resnorm, residual] = lsqnonneg(k2, l2)
    dres = abs(resnorm - 0.3695)          # compare with matlab result
    print 'ok, diff:', dres
    if dres > 0.11:
        raise Exeption('Error')
 
#    [x, resnorm, residual] = lsqnonneg(k1, l)
#    dres = abs(resnorm - 2.8639)          # compare with matlab result
#    print 'ok, diff:', dres
#    if dres > 0.01:
#        raise Exeption('Error')
 
#    C = numpy.array([[1.0, 1.0],
#                     [2.0, 1.0],
#                     [5.0, 1.0],
#                     [6.0, 1.0],
#                     [10.0, 1.0]])
 
#    d = numpy.array([3, 5, 11, 13, 21])
# 
#    [x, resnorm, residual] = lsqnonneg(C, d)

    print [x, resnorm, residual]
