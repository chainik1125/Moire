import numpy as np
import makeH
from variables import *
import json


def convert_to_json(data_dic):
    json_dic={}
    for i in data_dic.keys():
        json_dic[str(i)]=data_dic[i]
    return json_dic

#Parameters
theta=1.05
phi=2*np.pi/3
chosen_basis=makeH.generate_basis(number_of_particles=2,particle_dofno_dic=testpdd)
kinetic_diag=[[-v,{0:makeH.p0,1:makeH.px,2:makeH.p0}],[v,{0:makeH.p0,1:makeH.py,2:makeH.p0}]]
test_tun_dic=[[1,(0,1),{0:makeH.p0,1:makeH.t1,2:makeH.p0}],[1,(0,2),{0:makeH.p0,1:makeH.t2,2:makeH.p0}],[1,(0,3),{0:makeH.p0,1:makeH.t3,2:makeH.p0}]]
import sys

#parallel parameters
parts=int(sys.argv[2])
part=int(sys.argv[1])


def gen_Hk(kx,ky):
    return makeH.generate_H2(kx=kx,ky=ky,theta=theta,phi=phi,state_list=chosen_basis,diag_terms=kinetic_diag,tun_terms=test_tun_dic,UHK=UHK,mu=mu)



# print(np.linalg.eigh(hamiltonian(K2[0],K2[1]))[0])

# exit()


Nk1=100
# Define 2D chern number
# mu is chemical potential that determines filling, Nk is how many k points in each direction to integrate over
# ham_k is my own specific Hamiltonian class, replace it in the code where I marked with a function that gives you your H(kx, ky)
def chern(ham_k, Nk=Nk1, mu=0):
    kgrid={}
	# REPLACE BELOW ######
	#system = Ham3DCube(ham_k, Nx=0, Ny=0, Nz=1)
    ham = gen_Hk
	# IMPORTANT: my ham is a function of (kx, ky) that spits out a matrix for each (kx, ky)
	# just replace these lines with something that gives you this function (takes kx, ky, spits out matrix) for you
	######################
    a1=np.array([1/2,np.sqrt(3)/2])
    a2=np.array([-1/2,np.sqrt(3)/2])
    G1=2*np.pi*np.array([1,np.sqrt(3)/3])
    G2=2*np.pi*np.array([-1,np.sqrt(3)/3])
    # G1=2*np.pi*np.array([1,0])
    # G2=2*np.pi*np.array([0,1])
    # print(np.dot(G2,a2)/(2*np.pi))
    # exit()



    dk = 2*np.pi/Nk
    dk1=np.sqrt(np.dot(G1,G1))/Nk
    dk2=np.sqrt(np.dot(G2,G2))/Nk


    momenta=np.linspace(-1.5*np.pi,1.5*np.pi-dk,Nk)

    momenta1 = np.linspace(-QTR,QTR-QTR/Nk , Nk)
    momenta2 = np.linspace(-QTL, QTL-QTL/Nk, Nk)
    # momenta1 = np.linspace(0,G1-G1/Nk, Nk)
    # momenta2 = np.linspace(0, G2-G2/Nk, Nk)
    dof=len(gen_Hk(0,0))#Already defined in the generateHamilotnian file

    #dof = ham_k.dof() # this is the number of dof per unit cell (e.g. the bloch hamiltonian is dof x dof)
    states = dof

	#P is projector of occupied states, Q is unoccupied, dxP is derivative d_(kx) P
	#Layer resolved mag gives magnetization of each layer, sum over half to get surface mag
    P = np.zeros((states, states, Nk, Nk), dtype=complex)
    Q = np.zeros((states, states, Nk, Nk), dtype=complex)
    dxP = np.zeros((states, states, Nk, Nk), dtype=complex)
    dyP = np.zeros((states, states, Nk, Nk), dtype=complex)
	
	# Chern number
    C = 0

	# run through all k-space points to create projector of filled states 
    count=0
    # filled = (all states with energy below chemical potential mu)
    part_Nk=np.array_split(np.asarray(range(Nk)),parts)[part-1]#Just need to be careful that you don't get parts copies of the boundary contrib.

    #print(np.array_split(np.asarray(range(Nk)),parts))
    percent=int(0)
    for x in part_Nk:
        if int(len(part_Nk)/100)!=percent:
            percent=int(len(part_Nk))
            print(f"{percent}% done")
        # if x%(len(part_Nk)/100)==0:
        #     count+=1
        #     print(f"Projector constructing {count}% done")
        for y in range(Nk):
			
            k1 = momenta1[x][0]+momenta2[y][0]
            k2 = momenta1[x][1]+momenta2[y][1]
            # k1=momenta[x]
            # k2=momenta[y]
		
            vals, vecs = np.linalg.eigh(ham(k1, k2))
            vals = np.array(vals.real, dtype=float)
            args = vals.argsort()
            vals, vecs = vals[args], vecs[:,args]

            vals_filled = [vals[4]]#[e for e in vals if np.isclose(e,vals[0])]#[vals[0]]#Just pick the lowest energy state [vals[0]]
            # if len(vals_filled)>1:
            #     print(k1g,k2g)
            filled = len(vals_filled)
            vecs_filled = vecs[:, :filled]#Is vecs filled just a 
            P[:,:, x, y] = vecs_filled[:,:] @ vecs_filled[:,:].conj().T
            Q[:,:, x, y] = np.eye(states) - P[:,:, x, y]
	#run through k-space again to create derivatives of projection operators
    count=0
    for x in part_Nk:
        if x%(len(part_Nk)/100)==0:
            count+=1
            print(f"Derivative constructing {count}% done")
        for y in range(Nk):
            k1 = momenta1[x][0]+momenta2[y][0]
            k2 = momenta1[x][1]+momenta2[y][1]
            # k1=momenta[x]
            # k2=momenta[y]
		
			#I handled the boundary cases automatically by doing addition/subtraction (mod Nk)

            dxP[:,:,x,y] = (P[:,:, (x+1)%Nk, y] 
                - P[:,:, (x-1)%Nk, y]) / (2*dk)

            dyP[:,:,x,y] = (P[:,:, x, (y+1)%Nk] 
                - P[:,:, x, (y-1)%Nk]) / (2*dk)
                
            # if part==parts:
            #     dxP[:,:,x,y] = (P[:,:, (x+1), y] 
            #         - P[:,:, (x-1), y]) / (2*dk1)
            #     dyP[:,:,x,y] = (P[:,:, x, (y+1)] 
            #         - P[:,:, x, (y-1)]) / (2*dk2)
            # else:
            #     dxP[:,:,x,y] = (P[:,:, (x+1), y] 
            #         - P[:,:, (x-1), y]) / (2*dk1)
            #     dyP[:,:,x,y] = (P[:,:, x, (y+1)] 
            #         - P[:,:, x, (y-1)]) / (2*dk2)
                 
            # kgrid[(k1,k2)]=np.real(np.trace(dxP[:,:,x,y]))
            F = (Q[:,:, x,y] @ dyP[:,:, x,y] @ dxP[:,:, x,y] @ Q[:,:, x,y])
            
            

			# Make the integrand
            int_combo_chern = -2*F * 2*np.pi

            # add the integrand

            D = -1/(2*np.pi)**2 * dk**2 * np.imag(np.trace(int_combo_chern))#Ah - is there some Jacobian becaus you've got non-orthogonal vectors?
            C+=D
            kgrid[(k1,k2)]=D
            
            

	
    return C,kgrid

testC,testG=chern(ham_k=None)
print(testC)
file_name=f"Cherngridnobzpart{part}of{parts}Nk{Nk1}"
with open(f"./{file_name}.json", "w+") as f:
        json.dump(convert_to_json(testG), f)


