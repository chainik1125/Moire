import makeH
import pathdiag
import numpy as np
import sys


#Global variables
from variables import *

start=str(sys.argv[1])
end=str(sys.argv[2])
    
ordered_basis1=makeH.generate_basis(number_of_particles=1,particle_dofno_dic=testpdd)
ordered_basis2=makeH.generate_basis(number_of_particles=2,particle_dofno_dic=testpdd)
ordered_basis3=makeH.generate_basis(number_of_particles=3,particle_dofno_dic=testpdd)
ordered_basis4=makeH.generate_basis(number_of_particles=4,particle_dofno_dic=testpdd)
ordered_basis8=makeH.generate_basis(number_of_particles=8,particle_dofno_dic=testpdd)



tqs=[np.array([-1,-1]),np.array([1,0]),np.array([0,1])]
testnonlayer={'sublattice':2,'spin':2}
shell_basis=makeH.generate_shell_basis(shell_count=2,q_vecs=tqs,number_of_particles=1,nonlayer=testnonlayer)
print(shell_basis[0].particle_dic[1].dof_dic)
print(makeH.eyepreservation(shell_basis[0]))
print(len(shell_basis))



kinetic_diag=[[-v,{0:makeH.p0,1:makeH.px,2:makeH.p0}],[v,{0:makeH.p0,1:makeH.py,2:makeH.p0}]]
test_tun_dic=[[1,(0,1),{0:makeH.p0,1:makeH.t1,2:makeH.p0}],[1,(0,2),{0:makeH.p0,1:makeH.t2,2:makeH.p0}],[1,(0,3),{0:makeH.p0,1:makeH.t3,2:makeH.p0}]]

def gen_Hk(kx,ky):
    return makeH.generate_H2(kx=kx,ky=ky,theta=theta,phi=phi,state_list=shell_basis,diag_terms=kinetic_diag,tun_terms=test_tun_dic,UHK=UHK,mu=mu)


UHK=0
mu=0
Utau=0
Uff=0
start='KM'
end ='GammaM'
pathdiag.chained_path_plot_link(path_list=[vars_dic[start],vars_dic[end]],kpoints=100,generate_Hk=gen_Hk,UHK=UHK,mu=mu,Utau=Utau,Umu=Umu,Uff=Uff,names_reversed=names_reversed)





