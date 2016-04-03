# Reference to this code: http://arxiv.org/pdf/1509.04298.pdf
# Author: Nicola Pancotti @ Max Planck Institute of Quantum Optics
# Email:  nicola.pancotti@mpq.mpg.de
# Date:   14 Sep 2015

# One can run this code by simply typing 
# python toffoli.pub.py 25 .5 .5 
# the three arguments are:
# 1: learning rate step
# 2-3: weights of the initial state of the ancillas
#      |A> = $2 * |1> + $3 * |2>
# In principle, one could also optimize over these parameters 

# Please note, that the choice of the 'learning rates'
# is determinant for the convergence of the algorithm. 
# In this code the learning rate is named after 'step';
# it is the first argument passed at the command line.

# The authors found useful to run many instances in parallel 
# of the code with different learning rates in a cluster.
# They also noticed that the magnitude of the learning rates
# decreases exponentially with the number of qubits.

import sys
import os

#The following line is useful only if you can NOT install the 
#libraries globally in the machine you are using (e.g. cluster) 

#sys.path.append("path/to/the/libraries")

import qutip as qt
import scipy as scp
import numpy as np
import functions as fp
from math import log, sqrt, sin, cos, pi
from cmath import exp
from random import randrange
from datetime import datetime


####################################################
#                GATE and DIMENSION
####################################################
N = 4        # n of qubits

T = qt.toffoli() # Gate to simulate: toffoli
print("I am simulating the Toffoli gate")

#F = qt.fredkin() # fredkin 
#print("I am simulating the Fredkin gate")

try:
    H = fp.Htof
    G = T
except NameError:
    H = fp.Hfred
    G = F

####################################################
#                USEFUL DEFINITIONS
####################################################

dimG = G.shape[0]  
eta = float(sys.argv[2]) # We assume all the initial states of the ancilla  
xi = float(sys.argv[3]) # the same and separable. One can insert the coefficient 
                         # from command line.

dCS = cos(eta)*qt.basis(2,0) + exp(1j*xi)*sin(eta)*qt.basis(2,1) #Initial state of the ancillas


MAX_NUM_STEPS=500 #Max number of optimization steps, eventually included by command line
#MAX_NUM_STEPS=int(sys.argv[4])

###########################################
#       Stochastic Gradient Descent
###########################################

step = float(sys.argv[1]) #Learning rate


# OUTPUT file for the data
outputdir='output'
if not os.path.exists(outputdir):
    os.makedirs(outputdir)
outfile = open(str(outputdir)+'/toffoli_'+str(step)+'coef_'+str(eta)+str(xi), 'w+')

#Initial random guess for the coefficients of the Hamiltoninan 
J = np.random.rand(13)

######################################################
#   Initial guess optimal from the paper
#####################################################
###### Optimal parameters for Fredkin gate #####
#J = np.array([-4.712, 13.60, 11.15, 8.400, pi, 1.025])

###### Optimal parameters for toffoli gate #####
'''
J=np.zeros(13)
J[6] = 15.1043
J[5] = -5.66207
J[1] = -8.95113
J[3] = -4.9753
J[9] = -19.0948
J[11] = -4.30184
J[8] = -2.43
J[10] = -4.98238
J[12] =-0.208699
'''

#############################################################


# Check the fidelity function before optimization:
# quality of the initial state
#Ham = H(J,N)
#F = fp.Fidelity(Ham,G,N,dCS)
#print("The initial Fidelity is: "+str(F))


# Opitimization parameters 
delta = 0.001 # : for the derivative 
check = 0      # : to check convergence
s = 0          # : intermediate value of the optimization 
time = step    # : time steps

startTime = datetime.now() # : to check running time
while True :
    
    time += step
    
    rho_0 = qt.rand_ket(N = 8, dims = [[2,2,2], [1,1,1]]) # Initial random state of the system
    index = randrange(len(J)) # pick the derivative direction: We do not compute the gradient.
                              # Rather, each step, we compute 
                              # the derivative with respect to a random direction 
    
    JdJ = [x for x in J]
    JdJ[index] += delta
    
    HJ = H(J,N) # We compute the Hamiltonian for J 
    HJdJ = H(JdJ,N) #and for a little perturbation of J: JdJ
    
    grad = (fp.Likelihood3(JdJ, rho_0,G,dCS,HJdJ) - fp.Likelihood3(J, rho_0,G,dCS,HJ))/delta # The gradient of the Likelihood
    J[index] = J[index] + grad/sqrt(time) # Updating step 
    
    # One could drop the following 4 lines if not interested in the optimization behaviour 
    HJ = H(J,N)
    s = fp.Likelihood3(J,rho_0,G,dCS,HJ) 
    outfile.write(str(time/step)+'    '+str(s)+'\n')        
    #print(str(time/step)+'    '+str(s))        
    #print(str(time/step)+'    '+str(fp.Fidelity(HJ,G,N,dCS)))        
    
    # Check for convergence 
    if abs(s) > 0.99 :
        check += 1    
        if check == 20:
            break
                        
    # Check for max number of steps reached        
    if time/step > MAX_NUM_STEPS :
        break
            
    Jopt = J
            
outfile.write('#'+str(Jopt)+'\n')        
outfile.write('#'+str(datetime.now()-startTime)+'\n')        
                
#print(datetime.now()-startTime)
#print(Jopt)
outfile.close()

