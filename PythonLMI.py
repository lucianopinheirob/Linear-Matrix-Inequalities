"""
Created on Tue Sep 14 11:06:35 2021

@author: Luciano Pinheiro Batista 173096
Questão 1.11 item d) Lista 1 IA892
"""
import cvxpy as cp 
import numpy as np
import time as time
import scipy as sp

"Formulação do problema"
Q = np.eye(4)

P = cp.Variable((4,4), symmetric=True)

A = cp.Parameter((4,4))
A.value = np.matrix('0.7678 0.0208 -0.4699 0.3017; -0.048 -0.0336 0.2512 0.4410; 0.0401 -0.2927 0.5710 -0.5605; 0.3751 0.1035 -0.2433 -0.5067')
 
F = cp.atoms.affine.bmat.bmat([[P-Q, A.T@P], [P@A, P]]) #Função utilizada para construir matriz em bloco
     
des1 = F
des2 = P

constraints = [des1 >> 0, des2 >> 0]

obj = cp.Minimize(cp.trace(P))

prob = cp.Problem(obj, constraints)

"Resolve a LMI e calcula o tempo de resolução"
tempo = time.time()
prob.solve(solver=cp.CVXOPT)
tempo = time.time()-tempo
print('P = \n', P.value)
print('\nTempo de resolução da LMI:', np.around(tempo,4), 'segundos')


"Calcula a norma do erro em relação a solução da equação de Lyapunov"
Plyap = sp.linalg.solve_discrete_lyapunov(np.transpose(A.value), Q, method=None)
norma_erro = np.linalg.norm(Plyap-P.value, ord=2)
print('\nNorma do erro em relação a solução da equação de Lyapunov:', np.around(norma_erro,8))

