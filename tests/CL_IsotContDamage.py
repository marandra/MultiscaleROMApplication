import numpy as np

#INPUTS·······································································
#       eps_n -------- Strain step n (Time t)
#       eps_n1 ------- Strain step n+1 (Time t+Dt)
#       r_n ---------- Strain-like internal variable step n (Time t)
#   Mat.Properties:
#       C_elastic ---- Constitutive Elastic Tensor
#       H ------------ Hardening moduli
#
#OUTPUTS
#       r_n1 --------- Updated(or not) value of r for step n+1 (Time t+Dt)
#       sigma_n1 ----- Cauchy stress vector step n+1 (Time step t+Dt)
#       C_tang_n1 ---- Tangent constitutive operator step n+1 (Time step t+Dt)
#··············································································


eps_n = 
eps_n1 =
r_n =
C_elastic =
H =
r_0 = (...)


def StrainNorm_OnlyT(sigma_eff, eps):    #Damage Criterion Model (Damage Surface)
    argument = eps @ sigma_eff
    return np.sqrt(argument)

def HardeningLaw(r):
    q_0 = r_0
    q_n1 = q_0 + H * (r - r_0) 
    return q


# 1 - Compute effective stress vector and strain norm
sigma_eff_n1 = C_elastic @ (eps_n1.transpose)
Tao_eps_n1   = StrainNorm_OnlyT(sigma_eff_n1, eps_n1)
       


# 2 - Itentify state at step n+1 (time t+Dt) and update (or not) variables

if Tao_eps_n1 <= r_n :  #ELASTIC (UNDAMAGING) PROCESS
    
    r_n1 = r_n      #The strain-like internal variable keeps its (maxhistorical) value
    q_n1 = HardeningLaw(r_n1) #The stress-like int var is also updated
    d_n1 = 1 - (q_n1/r_n1)               #The damage variable is also updated
    
    sigma_n1 = (1-d_n1)*sigma_eff_n1 
    C_tang_n1 = (1-d_n1)*C_elastic 


else:                   #(DAMAGING)LOADING PROCESS
    
    r_n1 = Tao_eps_n1   #The strain-like int var is uptdated to it's new maxhist value 
    q_n1 = HardeningLaw(r_n1) #The stress-like int var is also updated
    d_n1 = 1 - (q_n1/r_n1)               #The damage variable is also updated
    
    sigma_n1 = (1 - d_n1)*sigma_eff_n1   #Updated stresses
    C_tang_n1= (1 - d_n1)*C_elastic - ((q_n1-H*r_n1)/(r_n1^3))*(sigma_eff_n1@sigma_eff_n1.transpose)


























