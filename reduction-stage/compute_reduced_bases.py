import sys
import configparser
import glob
import numpy as np

#import scipy as sci


def lsqnonneg(C, d, x0=None, tol=None, itmax_factor=10):
    '''Linear least squares with nonnegativity constraints.
    (x, resnorm, residual) = lsqnonneg(C,d) returns the vector x that minimizes norm(d-C*x)
    subject to x >= 0, C and d must be real
    '''
 
    eps = 2.22e-16    # from matlab

    def norm1(x):
        return abs(x).sum().max()

 
    def msize(x, dim):
        s = x.shape
        if dim >= len(s):
            return 1
        else:
            return s[dim]


    if tol is None:
        tol = 10 * eps * norm1(C) * (max(C.shape) + 1)
    C = np.asarray(C)
    (m,n) = C.shape
    P = np.zeros(n)
    Z = np.arange(1, n + 1)
    if x0 is None:
        x = P
    else:
        if any(x0 < 0):
            x = P
        else:
            x = x0
    ZZ = Z
    resid = d - np.dot(C, x)
    w = np.dot(C.T, resid)
    outeriter = 0
    it = 0
    itmax = itmax_factor * n
    exitflag = 1

    # outer loop to put variables into set to hold positive coefficients
    while np.any(Z) and np.any(w[ZZ - 1] > tol):
        outeriter += 1
        t = w[ZZ - 1].argmax()
        t = ZZ[t]
        P[t - 1] = t
        Z[t - 1] = 0
        PP = np.where(P != 0)[0] + 1
        ZZ = np.where(Z != 0)[0] + 1
        CP = np.zeros(C.shape)
        CP[:, PP - 1] = C[:, PP - 1]
        CP[:, ZZ - 1] = np.zeros((m, msize(ZZ, 1)))
        z = np.dot(np.linalg.pinv(CP), d)
        z[ZZ - 1] = np.zeros((msize(ZZ, 1), msize(ZZ, 0)))

        # inner loop to remove elements from the positve set which no longer belong
        while np.any(z[PP-1] <= tol):
            it += 1
            if it > itmax:
                max_error = z[PP - 1].max()
                raise Exception('Exiting: Iteration count (=%d) exceeded\n Try raising the \
                                 tolerance tol. (max_error=%d)' % (it, max_error))
            QQ = np.where((z <= tol) & (P != 0))[0]
            alpha = min(x[QQ] / (x[QQ] - z[QQ]))
            x = x + alpha * (z - x)
            ij = np.where((abs(x) < tol) & (P != 0))[0] + 1
            Z[ij - 1] = ij
            P[ij - 1] = np.zeros(max(ij.shape))
            PP = np.where(P != 0)[0] + 1
            ZZ = np.where(Z != 0)[0] + 1
            CP[:, PP - 1] = C[:, PP - 1]
            CP[:, ZZ - 1] = np.zeros((m, msize(ZZ, 1)))
            z=np.dot(np.linalg.pinv(CP), d)
            z[ZZ - 1] = np.zeros((msize(ZZ, 1), msize(ZZ, 0)))
        x = z
        resid = d - np.dot(C, x)
        w = np.dot(C.T, resid)
    return (x, sum(resid * resid), resid)


def ComputeJandb(Modes, weights, factorLEQ=None):

    eps = 2.22e-16    # from matlab

    if factorLEQ is None:
        factorLEQ = 1.0
    # Exact integral - numerical integration
    INTexact = np.dot(Modes.T,weights)

    # Total microscale volume
    vol = np.sum(weights)
    sqrtVol = np.sqrt(weights)

    # Matrix of modified modes (with zero integral)
    Xf = np.zeros(Modes.shape)

    # Loops over the initial modes
    Xf = np.subtract(Modes, INTexact / vol)
    Xf = np.multiply(Xf.T, sqrtVol).T

    # Singular Value Decomposition
    [Lambda,SValues,VValues] = np.linalg.svd(Xf,full_matrices=False)

    # fixed tolerance to define the reduced modified set of modes
    tol = np.max(Modes.shape) * eps * np.max(SValues)
    RankXf = sum(i > tol for i in SValues)
    Lambda = Lambda[:,0:RankXf]
    J = Lambda.T
    Jw = factorLEQ * sqrtVol / np.sqrt(vol)

    # Adding last row related with the sqrt of gauss integration weigths
    J = np.vstack([J, Jw.T])

    # Initializing the RHS vector for the optimization problem
    b = np.append(np.zeros(Lambda.T.shape[0]), factorLEQ * np.sqrt(vol))
    return (J, b, INTexact)


def ComputeROQ(Modes, weights, factorLEQ, nGP, tol):
    # computation of integration points weights
    #[J, b, INTexact] = ComputeJandb(k, GaussWeights, factorLEQ)
    [J, b, INTexact] = ComputeJandb(Modes, weights, factorLEQ)
    M=len(weights)
    y = np.arange(M)

    # Resudual vector, initial guess
    r = b

    #sys.exit()

    # Number of iterations
    it = 0
    mPOS=0
    z = []
    Jnorm = np.sqrt(sum(np.multiply(J,J),0))

    #print(J.shape)

    #sys.exit()

    # Point Selection Algorithm
    while (np.linalg.norm(r) / np.linalg.norm(b) > tol and mPOS <= nGP):

        # 1. Compute new point
        ObjFun = np.dot((J[:, y]).T, r)
        div = np.multiply(Jnorm[y], np.linalg.norm(r))
        ObjFun = np.divide(ObjFun, div)
        s = ObjFun.argmax()
        t = y[s]

        # 2. Move i from set y to set z
        z = (np.append(z, t)).astype(int)
        y=np.delete(y, s)

        # 3. solving LS conventional problem
        x = np.linalg.lstsq(J[:, z], b)[0]
        if any(x < 0):
            # 3. Determime alpha for solving a NNLS
            [x, resnorm, residual] = lsqnonneg(J[:,z], b)

        # 3. Determime alpha for solving a NNLS
        #[x, resnorm, residual] = lsqnonneg(J[:,z], b)

        # 4. Update the residual
        r = b - np.dot(J[:,z],x)

        # 5. Update mPOS and k
        mPOS = len(np.where(x>0)[0])

        # Iteration counter
        it = it + 1
        print("k =", it, "--- mPOS =",mPOS, "--- error (%) = ", np.linalg.norm(r)/np.linalg.norm(b)*100)

    # 6. Postprocess of points - neglecting null weights
    INDzero = np.where(x == 0)[0]
    if any(INDzero):
        z = np.delete(z, INDzero)
    w = np.multiply(x, np.sqrt(weights[z]))

    print("Reduced Weights")
    print(w)
    #
    print("sum of reduced weights")
    print(np.sum(w))
    #
    print("GP's index")
    print(z)

    return(w, z)


#######################################
# Main
#######################################
if __name__=='__main__':

   # get parameters
   fname = sys.argv[1]
   conf = configparser.ConfigParser()
   conf.read(fname)
   file_format = conf['Parameters']['file_format']
   nr_elastic_snapshots = int(conf['Parameters']['nr_elastic_snapshots'])
   nr_elements = int(conf['Parameters']['nr_elements'])
   nr_integration_points = int(conf['Parameters']['nr_integration_points'])
   nr_strain_components = int(conf['Parameters']['nr_strain_components'])
   nr_energy_reduced_modes = int(conf['Parameters']['nr_energy_reduced_modes'])
   energy_file_name = conf['Parameters']['energy_file_name']
   strain_file_name = conf['Parameters']['strain_file_name']
   gaussweights_file_name = conf['Parameters']['gaussweights_file_name']
   inelasticity_flag_name = conf['Parameters']['inelflag_file_name']

   # get files
   trajectory_paths = glob.glob('*_?')
   energy_elastic_files = []
   strain_elastic_files = []
   energy_inelast_files = []
   strain_inelast_files = []

   gaussweights_file    = glob.glob(trajectory_paths[0] + '/' + gaussweights_file_name)

   # TODO : define the number of snapshots per trajectory
   inelastic_flag_global = np.empty([len(trajectory_paths), 61])

   iTest=0
   inelastic_flag_max = np.empty([len(trajectory_paths),1])
   for path in trajectory_paths:
       inelasticity_flag_files = []
       inelasticity_flag_files.extend(sorted(glob.glob(path + '/' + inelasticity_flag_name + '*')))
       for i, file in enumerate(inelasticity_flag_files):
           inelastic_flag_global[iTest, i] = np.fromfile(file, dtype=np.float32)
       inelastic_flag_max[iTest]=inelastic_flag_global[iTest, :].argmax()
       iTest=iTest+1

   #print(inelastic_flag_global)
   #print(inelastic_flag_max)

   # TODO : for the future, take into account loadin/unloading trajectory cases
   #for path in trajectory_paths:
   for j, path in enumerate(trajectory_paths):
       #energy_elastic_files.extend(sorted(glob.glob(path + '/' + energy_file_name + '*'))[:nr_elastic_snapshots])
       #strain_elastic_files.extend(sorted(glob.glob(path + '/' + strain_file_name + '*'))[:nr_elastic_snapshots])
       #energy_inelast_files.extend(sorted(glob.glob(path + '/' + energy_file_name + '*'))[nr_elastic_snapshots:])
       #strain_inelast_files.extend(sorted(glob.glob(path + '/' + strain_file_name + '*'))[nr_elastic_snapshots:])
       #print(trajectory_paths[iPath])
       nr_elastic_snapshots=int(inelastic_flag_max[j][0])
       energy_elastic_files.extend(sorted(glob.glob(path + '/' + energy_file_name + '*'))[:nr_elastic_snapshots])
       strain_elastic_files.extend(sorted(glob.glob(path + '/' + strain_file_name + '*'))[:nr_elastic_snapshots])
       energy_inelast_files.extend(sorted(glob.glob(path + '/' + energy_file_name + '*'))[nr_elastic_snapshots:])
       strain_inelast_files.extend(sorted(glob.glob(path + '/' + strain_file_name + '*'))[nr_elastic_snapshots:])
       #print(strain_elastic_files)
       #print(strain_elastic_files)
       #print(strain_inelast_files)

   #print(trajectory_paths)

   #for ener in inelasticity_flag_files:
   #    print(ener)

   print("****************")
   print("STRAIN SNAPSHOTS")
   print("****************")
   # TODO: incluir el SVD del paquete scipy, parece que es mas optimo, ver si esta instalada una version reciente de scipy en el cluster
   # first part: read and compute elastic strain modes
   nr_dofs = nr_elements * nr_integration_points * nr_strain_components
   X = np.empty([nr_dofs, len(strain_elastic_files)])
   for i, file in enumerate(strain_elastic_files):
       #print(file)
       X[:, i] = np.fromfile(file, dtype=np.float32)
   [Utmp,Stmp,Vtmp] = np.linalg.svd(X, full_matrices=False)#[0]
   #print(Stmp)
   print("Step 01: svd of strain elastic snapshots")

   tolSVD=1e-8
   #tolSVD=1e-6 # for 27k 3D microscale
   cont=1
   Us_el=[]
   for iValue, sin_val in enumerate(Stmp):
       #print(Stmp[iValue])
       if sin_val>tolSVD:
           if cont==1:
               Us_el=Utmp[:,iValue]
               #print(Us_el)
           else:
               Us_el = np.column_stack((Us_el,Utmp[:,iValue]))
               #print(Us_el[:,iValue])
           cont=cont+1
   #print(Us_el.shape)
   print("Step 02: selection process of strain elastic modes")

   #[Us_el,ss_el,Vs_el] = np.linalg.svd(X, full_matrices=False)#[0]
   #Us_el = np.linalg.svd(X, full_matrices=False)[0]
   #[Us_el,ss_el,Vs_el] = sci.linalg.svd(X, full_matrices=False)#[0]
   #print(ss_el)
   #print(Us_el[:,0])
   #print(Us_el)

   # second part: read inelastic strain modes, remove elastic component, decomp svd
   X = np.empty([nr_dofs, len(strain_inelast_files)])

   #for i, file in enumerate(strain_inelast_files):
   #    X[:, i] = np.fromfile(file, dtype=np.float32)
   #print("Antes de proyector")
   #X = X - np.dot(Us_el, np.dot(Us_el.T, X))


   for i, file in enumerate(strain_inelast_files):
       #print(i)
       X[:, i] = np.fromfile(file, dtype=np.float32)
       for j in range(Us_el.shape[1]):
           X[:,i] = X[:,i] - np.multiply(np.dot(Us_el[:,j],X[:,i]),Us_el[:,j])

   #print("Desp de proyector")
   print("Step 03: projection of strain inelastic snapshots")
   #Us_in = np.linalg.svd(X, full_matrices=False)#[0]
   [Utmp,Stmp,Vtmp] = np.linalg.svd(X, full_matrices=False)#[0]
   #print(Stmp)
   print("Step 04: svd of strain inelastic modified snapshots")


   tolSVD=1e-7
   cont=1
   Us_in=[]
   #for iValue in range(len(Stmp)):
   for iValue, sin_val in enumerate(Stmp):
       #print(Stmp[iValue])
       if sin_val>tolSVD:
           if cont==1:
               Us_in=Utmp[:,iValue]
           else:
               Us_in = np.column_stack((Us_in,Utmp[:,iValue]))
           cont=cont+1
   #print(Us_in.shape)
   print("Step 05: Selection process of strain inelastic modes")

   Us=np.hstack([Us_el, Us_in])
   print("Step 06: assembly of global matrix of strain modes")
   print("SIZE OF Us")
   print(Us.shape)

   #sys.exit()

   print("****************")
   print("ENERGY SNAPSHOTS")
   print("****************")
   # third part: read and compute elastic energy modes, compute projector
   nr_dofs = nr_elements * nr_integration_points
   X = np.empty([nr_dofs, len(energy_elastic_files)])
   for i, file in enumerate(energy_elastic_files):
       X[:, i] = np.fromfile(file, dtype=np.float32)
   #Ue_el = np.linalg.svd(X, full_matrices=False)[0]
   [Utmp,Stmp,Vtmp] = np.linalg.svd(X, full_matrices=False)#[0]
   #print(Stmp)
   #print(Ue_el.shape)

   print("Step 07: svd of elastic energy snapshots")

   tolSVD=1e-7
   #tolSVD=1e-6 #for 27k 3D microscale
   cont=1
   Ue_el=[]
   for iValue, sin_val in enumerate(Stmp):
       #print(Stmp[iValue])
       if sin_val>tolSVD:
           if cont==1:
               Ue_el=Utmp[:,iValue]
           else:
               Ue_el = np.column_stack((Ue_el,Utmp[:,iValue]))
           cont=cont+1
   #print(Ue_el.shape)
   print("Step 08: Selection process of elastic energy modes")

   #sys.exit()

   # fourth part: read inelastic energy modes, remove elastic component, decomp svd
   X = np.empty([nr_dofs, len(energy_inelast_files)])
   #for i, file in enumerate(energy_inelast_files):
   #    X[:, i] = np.fromfile(file, dtype=np.float32)
   #print("Antes de proyector")
   #X = X - np.dot(Ue_el, np.dot(Ue_el.T, X))

   for i, file in enumerate(energy_inelast_files):
       #print(i)
       X[:, i] = np.fromfile(file, dtype=np.float32)
       for j in range(Ue_el.shape[1]):
           X[:,i] = X[:,i] - np.multiply(np.dot(Ue_el[:,j],X[:,i]),Ue_el[:,j])

   print("Step 09: projection of energy inelastic snapshots")
   #print("Desp de proyector")
   #print(X.shape)
   #Ue_in = np.linalg.svd(X, full_matrices=False)[0]
   [Utmp,Stmp,Vtmp] = np.linalg.svd(X, full_matrices=False)#[0]
   #print(Stmp)
   print("Step 10: svd of inelastic energy modified snapshots")


   tolSVD=1e-7
   cont=1
   Ue_in=[]
   for iValue, sin_val in enumerate(Stmp):
       #print(Stmp[iValue])
       if sin_val>tolSVD:
           if cont==1:
               Ue_in=Utmp[:,iValue]
           else:
               Ue_in = np.column_stack((Ue_in,Utmp[:,iValue]))
           cont=cont+1
   #print(Ue_in.shape)
   print("Step 11: selection process of inelastic energy modes")

   Ue=np.hstack([Ue_el, Ue_in])
   print("Step 12: assembly of global matrix of energy modes")
   print("SIZE OF Ue")
   print(Ue.shape)

   #sys.exit()

   # reading the gauss weigths for computing the ROQ
   nr_dofs = nr_elements * nr_integration_points
   for i, file in enumerate(gaussweights_file):
       GaussWeights = np.fromfile(file, dtype=np.float32)
   #print("Gauss Weights")
   #print(GaussWeights)

   #sys.exit()

   # TODO: change the print format of the matrix in order to avoid wrong tabulation because of the minus (-) sign.
   print("**********************")
   print("Printing data to files")
   print("**********************")
   file_Us='./Output/StrainBasis.dat'
   with open(file_Us,'wb') as ofileBasis:
       #np.savetxt(ofileBasis, Us, fmt='%+.13f')
       np.savetxt(ofileBasis, Us, fmt='%.13f')
       #for v in Us_el:
           #ofileBasis.write("f\n".format(v))
           #np.savetxt(ofileBasis, v, fmt='%.13f')
           #np.matrix.tofile(ofileBasis, v, fmt='%.13f')

   #num_ener_modes=20
   if nr_energy_reduced_modes > Ue.shape[1]:
       sys.exit("Error: number of energy modes greater than the total number of computed energy modes")
   else:
       num_ener_modes = nr_energy_reduced_modes

   Ue_red=Ue[:,0:num_ener_modes]

   file_Ue='./Output/EnergyBasis.dat'
   with open(file_Ue,'wb') as ofileBasis_ener:
       #np.savetxt(ofileBasis_ener, Ue, fmt='%+.13f')
       np.savetxt(ofileBasis_ener, Ue_red, fmt='%.13f')

   #cc=np.array([0, 9, 19, 8, 6])
   #d=np.where(cc==7)[0]
   #if not d:
   #    print("no hay na!!!")
   #else:
   #    print(d)

   #sys.exit()

   print("************************")
   print("REDUCED ORDER QUADRATURE")
   print("************************")
   # fifth part: COMPUTING REDUCED ORDER QUADRATURE (ROQ)
   factorLEQ=1.0
   tol = 1e-10
   nGP = num_ener_modes
   #nGP = 20
   #nGP=Ue.shape[1] #In case of use the same number of points as energy modes.
   #nGP = 120
   #nGP = 12
   [w,z] = ComputeROQ(Ue_red, GaussWeights, factorLEQ, nGP, tol)

   #sys.exit()

   # Print matrix with new weigths
   roq_weigths=np.empty([nr_elements, nr_integration_points])
   for i in range(nr_elements):
       for j in range(nr_integration_points):
           i_elem=nr_integration_points*(i)+j
           d=np.where(z==i_elem)[0]
           if not d:
               roq_weigths[i][j]=-1
           else:
               roq_weigths[i][j]=w[d[0]]

   file_roq_weights='./Output/weights.dat'
   with open(file_roq_weights,'wb') as ofile_roq_weights:
       np.savetxt(ofile_roq_weights, roq_weigths, fmt='%+.13f')