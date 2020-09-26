# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 10:20:46 2020

@author: bguerle
"""

## This code has been written by Burak Gurlek for QC Mentorship program selection Task 1
##!!! Please use cirq==0.8.2 version

# import packages
import matplotlib.pyplot as plt
import numpy as np
import cirq
import sympy
import scipy

# Define circuit parameters 
qubits = cirq.LineQubit.range(4) # create 4 qubits
num_reps=100 # number of measurement of QC circuit
maxiter=1000 # number of iteration in minimization
L=6 # number of circuit layer to be aplied

# Define goal qubit vector
phi=np.round(np.random.rand(1,len(qubits))) # create random 4 - qubit array 

# =============================================================================
#  Define functions to be used
# =============================================================================

def Uodd_layer(angles,gatetype):
  """Generator for Uodd(Theta,gatetype) layer
     requires qubit rotation angles for input, dimension  1x4
     and gatetype which can be RX,RY or RZ
  """
  # requires qubit rotation angles for input, dimension  1x4
  
  if gatetype=='RX': # create rotation gate with angles[ii]) on qubits[ii]
     for ii in range(len(qubits)):
         yield cirq.rx(angles[ii])(qubits[ii])
  elif gatetype=='RY':
     for ii in range(len(qubits)):
         yield cirq.ry(angles[ii])(qubits[ii])
  elif gatetype=='RZ':
     for ii in range(len(qubits)):
         yield cirq.rz(angles[ii])(qubits[ii])
  else:
      print('Uodd Layer can only include RX or RY or RZ gate ')

def Ueven_layer(angles,gatetype):
      """Generator for Ueven(Theta,gatetype) layer
         requires qubit rotation angles for input, dimension  1x4
         and gate type which can be RX,RY or RZ
      """
      if gatetype=='RX': # create rotation gate with angles[ii]) on qubits[ii]
         for ii in range(len(qubits)):
             yield cirq.rx(angles[ii])(qubits[ii])
      elif gatetype=='RY':
           for ii in range(len(qubits)):
               yield cirq.ry(angles[ii])(qubits[ii])
      elif gatetype=='RZ':
           for ii in range(len(qubits)):
               yield cirq.rz(angles[ii])(qubits[ii])
      else:
           print('Ueven Layer can only include RX or RY or RZ gate before double qubit gates ')             
          
      for ii in range(len(qubits)): # create CZ gate between qubits[ii] & qubits[jj]
            for jj in range(ii+1,len(qubits)):
                yield cirq.CZ(qubits[ii],qubits[jj])
          
def QCircuit(qubits,L,Uodd_gatetype,Ueven_gatetype):
    """Generator for QC Circuit
       requires qubits as input, number of layers in the circuit
       and gate type for even and odd U layers
       returns QCircuit (Quantum Circuit)
    """
    optmcirc=cirq.Circuit()
    # Creating angle symbol array for assignmet to rotation gates: Th01, first index for layer, 
    # second one for qubit in the layer
    angles=([[sympy.symbols('Th%d%d'%(i,j)) for j in range(len(qubits))] for i in range(L)])
    
    for ii in range(L): # define layer and append to circuit
        optmcirc.append(Uodd_layer(tuple(angles[ii]),Uodd_gatetype),strategy=cirq.InsertStrategy.NEW_THEN_INLINE)
        optmcirc.append(Ueven_layer(tuple(angles[ii]),Ueven_gatetype),strategy=cirq.InsertStrategy.NEW_THEN_INLINE)
    
    optmcirc.append(cirq.measure(*[row for row in qubits],key='m'))
    return optmcirc

def distance(ii_angles,noflayers,circuit):
    """Calculation of distance between resulting measurements and goal qubit
       requires qubit rotation angles for input, dimension  1xL*4, number of layers
       and quantum circuit
       returns distance between measured qubits from simulation and the goal qubit (phi)
    """
    ii_angles=np.remainder(ii_angles, 2*np.pi) # wrap angles to [0,2pi]
    
    # convert angle array to matrix such that each row corresponds to angles of a layer
    i_angles=[]
    [i_angles.append(ii_angles[len(qubits)*i:len(qubits)*i+len(qubits)]) for i in range(L)]
    
    # assign angles to parametrized 'Th00' via dict
    param_dict={}
    for i in range(noflayers):
        for j in range(len(qubits)):
            param_dict['Th%d%d' %(i,j)]=i_angles[i][j]
    # create parametric solver
    params=cirq.ParamResolver(param_dict)
    
    # simulate the circuit
    sim=cirq.Simulator()
    results=sim.run(circuit, param_resolver=params,repetitions=num_reps)
    
    # unwrap simulation results to 4xnum_reps array, such that each column represnt
    # another measurement result
    raw_meas=str(results)[2:].split(", ")
    meas_qubits=[]
    for ii in range(num_reps):
        meas_qubits.append([])
        for jj in range(len(raw_meas)):
            meas_qubits[ii].append(int(raw_meas[jj][ii])) # measured qubit matrix
    
    # calculate distance between measured qubits (meas_qubits) and goal qubit (phi)  
    dist=0
    for ii in range(len(meas_qubits)):
        dist=dist+np.linalg.norm(meas_qubits[ii]-phi)
    dist=dist/num_reps
    return dist

def run_Trial(ii_angles,noflayers):
    """run Qcircuit with optimized rotation angles
       requires angle array in rad, dimension 1xL*4 and number of layers
       return measured_qubits, distance with goal qubit and results object
       given input rotation angles
    """
    ii_angles=np.remainder(ii_angles, 2*np.pi) # wrap angles to [0,2pi]
    
    # convert angle array to matrix such that each row corresponds to angles of a layer
    i_angles=[]
    [i_angles.append(ii_angles[len(qubits)*i:len(qubits)*i+len(qubits)]) for i in range(L)]
    
    # assign angles to parametrized 'Th00' etc via dict
    param_dict={}
    for i in range(noflayers):
        for j in range(len(qubits)):
            param_dict['Th%d%d' %(i,j)]=i_angles[i][j]
    # create parametric solver
    params=cirq.ParamResolver(param_dict)
    
    # simulate the circuit
    sim=cirq.Simulator()
    results=sim.run(circuit, param_resolver=params,repetitions=num_reps)
    
    # unwrap simulation results to 4xnum_reps array, such that each column represnt
    # another measurement result
    raw_meas=str(results)[2:].split(", ")
    meas_qubits=[]
    for ii in range(num_reps):
        meas_qubits.append([])
        for jj in range(len(raw_meas)):
            meas_qubits[ii].append(int(raw_meas[jj][ii]))
    
    # calculate distance between measured qubits (meas_qubits) and goal qubit (phi)  
    dist=0
    for ii in range(len(meas_qubits)):
        dist=dist+np.linalg.norm(meas_qubits[ii]-phi)
    dist=dist/num_reps
    
    return meas_qubits, dist, results
 
def callbackF(Xi):
    ''' Callback function that prints the progress of optimization algorithm'''
    global Nfeval
    global noflayers
    global circuit
    print ('{0:4d}  {1:3.10f}' .format(Nfeval, distance(Xi,noflayers,circuit))) # print evoluation number with the calculated distance between qubits
    Nfeval += 1


# =============================================================================
#  End of function definitions
# =============================================================================




# Create array for minimized results
res=[] # for optimization results object
meas_optm=[] # for measured states with optimized angles
meas_dist=[] # for measured distance with optimized distance
meas_results=[] # for measured results with optimized distance
optm_angles=[] # optimized angles first for for L=1 etc.
# minimization of the circuit for different # of Layers
for ii in range(1,L+1):
    
    # Create quantum circuit
    circuit=QCircuit(qubits,ii,'RX','RZ')
    print(circuit)

    # Convert initialized 2D radom matrix to 1D array
    x0=np.random.uniform(0, 2*np.pi, size=len(qubits)*ii) # initial angle array
    Nfeval=0 # number of evoluations of the optimization algorithm
    noflayers=ii
    print ('for L={0:d}' .format(ii))
           
    res.append(scipy.optimize.minimize(distance,x0,args=(ii,circuit),method='Powell',callback=callbackF,options={'xtol':1e-12, 'ftol':1e-12 ,'maxiter': maxiter}))
    if res[ii-1].success==True:
        print ('Minimum found')
        # storing data
        store_Trial=run_Trial(res[ii-1].x,ii)
        meas_optm.append(store_Trial[0]) # store measured states with optimized angles
        meas_dist.append(store_Trial[1]) # store measured distance with optimized angles
        meas_results.append(store_Trial[2]) # store measurement results with optimized angles
        optm_angles.append(np.remainder(res[ii-1].x, 2*np.pi)) # optimized angles wrapped in [0,2pi] 

    else:
        print('Convergence is not achieved, use different algorith or parameters')
# Plot most common 5 outputs from QC simulation with optimized andles
for ii in range(1,L+1):
    hist = meas_results[ii-1].histogram(key='m')
    num = 15
    probs = [v/meas_results[ii-1].repetitions for _,v in hist.most_common(num)]
    plt.figure(figsize=(9, 3))
    plt.title('Probability of {} Most Common Outputs for QC Circuit with L={} Layers'.format(num,ii))
    plt.bar([x for x in range(len(probs))],probs)
    plt.show()

# Plot distance between QC output and goal qubit vs number of layers in QC circuit
plt.figure(figsize=(9, 3))
plt.plot(range(1,L+1), meas_dist, 'bs')
plt.xlabel('No. of Layers in QC')
plt.ylabel('Distance ')
plt.show()

print('Sometimes with increase of L I do get divergent results, this happens when I do find the ground state with L=1. I guess there may be a problem with my optimizer since the dimension of the space increases, it could not find the right jacobian eaisly.It may stuck at local minima')
print('I would be glad if you can let me know the error I did')
# =============================================================================
# Check what other parametrized gates do e.g. RY for U1 and RX for U2
# =============================================================================

# =============================================================================
# 
# 
# # Create array for minimized results
# res_YX=[] # for optimization results object
# meas_optm_YX=[] # for measured states with optimized angles
# meas_dist_YX=[] # for measured distance with optimized distance
# meas_results_YX=[] # for measured results with optimized distance
# optm_angles_YX=[] # optimized angles first for for L=1 etc.
# # minimization of the circuit for different # of Layers
# for ii in range(1,L+1):
#     
#     # Create quantum circuit
#     circuit=QCircuit(qubits,ii,'RY','RX')
#     print(circuit)
# 
#     # Convert initialized 2D radom matrix to 1D array
#     x0=np.random.uniform(0, 2*np.pi, size=len(qubits)*ii) # initial angle array
#     Nfeval=0 # number of evoluations of the optimization algorithm
#     noflayers=ii
#     print ('for L={0:d}' .format(ii))
#            
#     res_YX.append(scipy.optimize.minimize(distance,x0,args=(ii,circuit),method='Powell',callback=callbackF,options={'xtol':1e-10, 'ftol':1e-12 ,'maxiter': maxiter}))
#     if res_YX[ii-1].success==True:
#         print ('Minimum found')
#         # storing data
#         store_Trial_YX=run_Trial(res[ii-1].x,ii)
#         meas_optm_YX.append(store_Trial[0]) # store measured states with optimized angles
#         meas_dist_YX.append(store_Trial[1]) # store measured distance with optimized angles
#         meas_results_YX.append(store_Trial[2]) # store measurement results with optimized angles
#         optm_angles_YX.append(np.remainder(res[ii-1].x, 2*np.pi)) # optimized angles wrapped in [0,2pi] 
# 
#     else:
#         print('Convergence is not achieved, use different algorith or parameters')
# # Plot most common 5 outputs from QC simulation with optimized andles
# for ii in range(1,L+1):
#     hist = meas_results_YX[ii-1].histogram(key='m')
#     num = 15
#     probs = [v/meas_results_YX[ii-1].repetitions for _,v in hist.most_common(num)]
#     plt.figure(figsize=(9, 3))
#     plt.title('Probability of {} Most Common Outputs for QC Circuit with L={} Layers [U1 Layer RY,U1 Layer RX ]'.format(num,ii))
#     plt.bar([x for x in range(len(probs))],probs)
#     plt.show()
# 
# # Plot distance between QC output and goal qubit vs number of layers in QC circuit
# plt.figure(figsize=(9, 3))
# plt.plot(range(1,L+1), meas_dist_YX, 'bs')
# plt.title('Distance vs No. of Layers [U1 Layer RY,U1 Layer RX ]'.format(num,ii))
# plt.xlabel('No. of Layers in QC')
# plt.ylabel('Distance ')
# plt.show()
# 
# 
# =============================================================================


