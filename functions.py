# Reference to this code: http://arxiv.org/pdf/1509.04298.pdf
# Author: Nicola Pancotti @ Max Planck Institute of Quantum Optics
# Email:  nicola.pancotti@mpq.mpg.de
# Date:   14 Sep 2015

import qutip as qt
from math import log


def Likelihood2(J, rho_0, G, dCS, H) : #Likelihood function for 2 qubits gates
    
    Ak = G*rho_0
    A = Ak*Ak.dag()
    rho = qt.tensor(rho_0, dCS)
    U = (-1j*H).expm()
    Btemp = U*rho
    Bk = Btemp*Btemp.dag()
    B = Bk.ptrace([0,1])
    out = (A*B).tr()

    return abs(out)

def Likelihood3(J, rho_0, G, dCS, H) : #Likelihood function for 2 qubits gates
    
    
    Ak = G*rho_0
    A = Ak*Ak.dag()
    rho = qt.tensor(rho_0, dCS)
    U = (-1j*H).expm()
    Btemp = U*rho
    Bk = Btemp*Btemp.dag()
    B = Bk.ptrace([0,1,2])
    out = (A*B).tr()

    return abs(out)



def Hfred(x,N) : # A possible Hamiltonian for the Fredkin gate
    k = 0
    H = 0
    sx = qt.sigmax()/2
    sy = qt.sigmay()/2
    sz = qt.sigmaz()/2
    Id = qt.qeye(2)

    for q in [sx, sy, sz] :
        temp = 0
        OpChain = [Id]*N
        OpChain[2] = q
        OpChain[1] = q
        temp += x[k]*qt.tensor(OpChain)
        H += temp 
    k+=1        

    for q in [sx,sz]:
        
        temp = 0
        OpChain = [Id]*N
        OpChain[2] = q
        OpChain[0] = q
        temp += x[k]*qt.tensor(OpChain)
    
        OpChain = [Id]*N
        OpChain[1] = q
        OpChain[0] = q
        temp += x[k]*qt.tensor(OpChain)
        k += 1
        H += temp 
    
    for q in [1,2]:
        temp = 0
        OpChain = [Id]*N
        OpChain[q] = sx
        OpChain[3] = sx
        temp += x[k]*qt.tensor(OpChain)
        H += temp 
    k+=1        

    temp = 0
    OpChain = [Id]*N
    OpChain[0] = sz
    temp += x[k]*qt.tensor(OpChain)
    H += temp 
    k += 1
    
    temp = 0
    OpChain = [Id]*N
    OpChain[3] = sx
    temp += x[k]*qt.tensor(OpChain)#last one

    H += temp 

    
    return H



def Htof(x,N) :
    
    k = 0
    H = 0

    sx = qt.sigmax()/2
    sz = qt.sigmaz()/2
    Id = qt.qeye(2)


    for q in [sx,sz]:    
        temp = 0
        OpChain = [Id]*N
        OpChain[0] = q
        OpChain[1] = q
        temp += x[k]*qt.tensor(OpChain)
        k+=1        
        H += temp 

    for p in [2,3]:
        for q in [sx,sz]:
            
            temp = 0
            OpChain = [Id]*N
            OpChain[0] = q
            OpChain[p] = q
            
            temp += x[k]*qt.tensor(OpChain)
            
            OpChain = [Id]*N
            OpChain[1] = q
            OpChain[p] = q
            
            temp += x[k]*qt.tensor(OpChain)
            k += 1
            
            H += temp 
            
   
    for q in [sx,sz]:
        
            
        temp = 0
        OpChain = [Id]*N
        OpChain[2] = q
        OpChain[3] = q
        temp += x[k]*qt.tensor(OpChain)
        k+=1    
        H += temp 
            
            
    for i in range(2) :
        
        temp = 0
        OpChain = [Id]*N
        OpChain[i] = sz
        temp += x[k]*qt.tensor(OpChain)
        H += temp 

    k += 1

    for i in [2,3] :
    
        temp = 0
        OpChain = [Id]*N
        OpChain[i] = sx
        temp += x[k]*qt.tensor(OpChain)
        k += 1    

        OpChain = [Id]*N
        OpChain[i] = sz
        temp += x[k]*qt.tensor(OpChain)#last one
        k += 1

        H += temp 

    
    return H


###############################################
#     FUNCTIONS TO IMPLEMENT THE FIDELITY
#     useful to double check the stochastic
#     optimization
##############################################



def getGate (G): #extract non zero elements of the gate and save them in s
    
    s = []

    rows = G.shape[0]
    colums = G.shape[1]
    
    for i in range(rows):
        for j in range(colums):
            if G[i][0][j] != 0 :
                s.append([i,j])
    return s

def binary(a):  #binary representation of a number a: useful to write the computational basis 
    s=''                                          #inside the Fidelity
    t={'0':'000','1':'001','2':'010','3':'011',
       '4':'100','5':'101','6':'110','7':'111'}
    for c in oct(a)[1:]:
            s+=t[c]
    return s

def getBasis (a,G) : #get the basis states according to the binary of a: 10 -> |10>
    
    dimG = G.shape[0]  
    CareStateDim = int(log(dimG,2))

    if a == 0:
          B = [qt.basis(2,0)]*(CareStateDim)
          return qt.tensor(B)
    
    c = binary(a)
    if dimG == 8:
        return qt.tensor(qt.basis(2,int(c[0])) , qt.basis(2,int(c[1])), qt.basis(2,int(c[2])))
    if dimG == 4:
        return qt.tensor(qt.basis(2,int(c[1])) , qt.basis(2,int(c[2])))
    
def Fidelity (Ham,G,N,dCS): #(J):  #Fidelity function
    
    s = getGate(G)
    #H = HamiltonianAB(J)
    dimG = G.shape[0]  
    Fid = 1./(dimG + 1)
    U = (-1j*Ham).expm()
    Udag = (1j*Ham).expm()
    
    for x in s :
        for y in s:
            
            #definition of the basis kets and bras.             
            bra_i = getBasis(x[0],G).dag()
            ket_j = getBasis(y[0],G)
            ket_k = getBasis(x[1],G)
            bra_l = getBasis(y[1],G).dag()
            
            Epsilon = U*qt.tensor(ket_k*bra_l, dCS*dCS.dag())*Udag
            Eps_ijkl = bra_i*(Epsilon.ptrace([0,1,2]))*ket_j
            
            Gstar_ik = G[x[0],x[1]].conjugate()
            G_jl = G[y[0],y[1]]
            
            
            fidStep = (1./(dimG*(dimG+1)))*Gstar_ik*Eps_ijkl*G_jl     
            Fid += fidStep[0][0][0]
            
    return abs(Fid)
