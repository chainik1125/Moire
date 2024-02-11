import numpy as np
import makeH
from variables import *
import itertools




#c6z=({0:makeH.p0,1:pexpz6},1)

c2x=({0:makeH.px,1:makeH.px},1)






# def tpp_symmetry_action(state_list,term):
#     #Do the same thing as you would for constructing the Hamiltonian matrix - the exception being that here you 
#     #act on all the particles at once, rather than individually
#     H0=np.zeros((len(state_list),len(state_list)),dtype=complex)
#     for state in state_list:
        
#     return None

def diagtpp(state_list,layer_pauli_dic,prefactor):

    H1=np.zeros((len(state_list),len(state_list)),dtype=complex)
    particle_dic_temp1={}#Idea is to constantly overwrite
    dof_dic_temp={}

    for state in state_list:
        result_list=[]
        for k in state.particle_dic.keys():
            coeff=1
            dof_dic_temp[0]=layer_pauli_dic[0](state.particle_dic[k].dof_dic[0])[1]
            coeff=coeff*layer_pauli_dic[0](state.particle_dic[k].dof_dic[0])[0]
            temp_k=kdic[state.particle_dic[k].dof_dic[0]]
            temp_thetak=get_theta(k=temp_k)-theta/2
            coeff_dic={px:np.cos(temp_thetak),py:np.sin(temp_thetak)}
            #sublattice and spin part
            for l in range(1,len(layer_pauli_dic.keys())):
                dof_dic_temp[l]=layer_pauli_dic[l](state.particle_dic[k].dof_dic[l])[1]
                coeff=coeff*layer_pauli_dic[l](state.particle_dic[k].dof_dic[l])[0]
            coeff=coeff_dic[layer_pauli_dic[1]]*np.sqrt(np.vdot(temp_k,temp_k))
            particlen=particle(dof_dic=dict(dof_dic_temp))
            # print(vars(particlen))
            particle_dic_temp1[k]=particlen
            for j in state.particle_dic.keys():
                if j!=k:
                    particle_dic_temp1[j]=state.particle_dic[j]
            new_state=basis(particle_dic=dict(particle_dic_temp1))
            # print(eyepreservation(new_state))
            result=(coeff,new_state)
            result_list.append(result)
        # for result in result_list:
                # print(f"Input state {eyepreservation(state)}")
                # print(f"Output state {eyepreservation(result[1])}")
                # print(len(result_list))
            # particle_dic_temp1.clear()
    #2.6 Now I need to map the states into some standard Hamiltonian ordering.
    #I need to do an unordered match of the state to a state in the basis, and then I need to count the number of swaps
#Have to do this because, since I've introduced the basis as a class - seperate instances of the class count as different objects...    
        for result in result_list:
            for i in state_list:
                if basis.eqnoorder(i,result[1]):
                    temp_H=np.zeros((len(state_list),len(state_list)),dtype=complex)
                    swapcount=basis.swaps(result[1],i)
                    position1=(state_list.index(i),state_list.index(state))
                    temp_H[position1]=((-1)**swapcount)*prefactor*result[0]
                    H1=H1+temp_H
        # position1=(ordered_basis.index(basis1),ordered_basis.index(state)) #in the matrix, row index is final state, col. index is final state
        #the SO coupling brings imaginary terms in, unfortunately...
        # bit of a hack but I'm just using the pair to seperate out the coefficent from the state label
            
        return H1
    
def tpp2(state_list,pauli_dic,prefactor):#for implementing symmetry transforms
    H1=np.zeros((len(state_list),len(state_list)),dtype=complex)
    #By convention this is the thing that will have the e^... in the exponent
    for state in state_list:
        coeff_list=[]
        dof_dic_temp1={}
        state_dic_temp={}
        for k in state.particle_dic.keys():
            coeff=1
            for l in state.particle_dic[k].dof_dic.keys():
                print(l)
                print(pauli_dic[l])
                print(state.particle_dic[k].dof_dic[l])
                print(makeH.eyepreservation(state))
                dof_dic_temp1[l]=pauli_dic[l](state.particle_dic[k].dof_dic[l])[1]
                coeff=coeff*pauli_dic[l](state.particle_dic[k].dof_dic[l])[0]    
            coeff_list.append(coeff)
            particlen=makeH.particle(dof_dic=dict(dof_dic_temp1))
            state_dic_temp[k]=particlen
        result_state=makeH.basis(particle_dic=dict(state_dic_temp))
            
        # particle4=particle(mu=pmu(state.p4.mu)[1],tau=ptau(state.p4.tau)[1],sigma=pz(state.p4.sigma)[1])
        # result4=pmu(state.p4.mu)[0]*ptau(state.p4.tau)[0]*np.exp(1j*psigmazexp*pz(state.p4.sigma)[0])
        # coeff_list.append(result4)
        for i in state_list:
            if makeH.basis.eqnoorder(result_state,i):
                temp_H=np.zeros((len(state_list),len(state_list)),dtype=complex)
                swapcount=makeH.basis.swaps(result_state,i)
                position1=(state_list.index(i),state_list.index(state))
                coeff1=1
                for k in coeff_list:
                    coeff1=coeff1*k
                temp_H[position1]=((-1)**swapcount)*prefactor*coeff1
                H1=H1+temp_H
    return H1

def tpp3(state_list,pauli_dic,prefactor):#for implementing symmetry transforms
    H1=np.zeros((len(state_list),len(state_list)),dtype=complex)
    #By convention this is the thing that will have the e^... in the exponent
    for state in state_list:
        coeff_list=[]
        dof_dic_temp1={}
        state_dic_temp={}
        for k in state.particle_dic.keys():
            coeff=1
            for l in state.particle_dic[k].dof_dic.keys():
                dof_dic_temp1[l]=pauli_dic[l](state.particle_dic[k].dof_dic[l])[1]
                coeff=coeff*pauli_dic[l](state.particle_dic[k].dof_dic[l])[0]    
            coeff_list.append(coeff)
            particlen=makeH.particle(dof_dic=dict(dof_dic_temp1))
            state_dic_temp[k]=particlen
        result_state=makeH.basis(particle_dic=dict(state_dic_temp))
            
        # particle4=particle(mu=pmu(state.p4.mu)[1],tau=ptau(state.p4.tau)[1],sigma=pz(state.p4.sigma)[1])
        # result4=pmu(state.p4.mu)[0]*ptau(state.p4.tau)[0]*np.exp(1j*psigmazexp*pz(state.p4.sigma)[0])
        # coeff_list.append(result4)
        for i in state_list:
            if makeH.basis.eqnoorder(result_state,i):
                temp_H=np.zeros((len(state_list),len(state_list)),dtype=complex)
                swapcount=makeH.basis.swaps(result_state,i)
                position1=(state_list.index(i),state_list.index(state))
                coeff1=1
                for k in coeff_list:
                    coeff1=coeff1*k
                temp_H[position1]=((-1)**swapcount)*prefactor*coeff1
                H1=H1+temp_H
    return H1


#monolayer stuff
# monopdd={'sublattice':2,'spin':2}
# mono_basis2=makeH.generate_basis(number_of_particles=2,particle_dofno_dic=monopdd)

# ordered_basis2=makeH.generate_basis(number_of_particles=2,particle_dofno_dic=makeH.testpdd)
# testind=3
# teststate=ordered_basis2[testind]
# testwf=np.zeros(len(mono_basis2),dtype=complex)
# testwf[testind]=1
# testsym=tpp2(state_list=mono_basis2,pauli_dic=c3z[0],prefactor=c3z[1])
# reswf=np.dot(testsym,testwf)
# for i in range(len(reswf)):
#     if np.isclose(np.abs(reswf[i]),0)==False:
#         print(f"Input state {makeH.eyepreservation(state=teststate)}, Output state {makeH.eyepreservation(state=ordered_basis2[i])}, amp={reswf[i]}")

# kinetic_diag=[[-v,{0:makeH.p0,1:makeH.px,2:makeH.p0}],[v,{0:makeH.p0,1:makeH.py,2:makeH.p0}]]
# test_tun_dic=[[1,(0,1),{0:makeH.p0,1:makeH.t1,2:makeH.p0}],[1,(0,2),{0:makeH.p0,1:makeH.t2,2:makeH.p0}],[1,(0,3),{0:makeH.p0,1:makeH.t3,2:makeH.p0}]]

# testH=makeH.generate_H2(kx=K[0],ky=K[1],theta=theta,phi=phi,state_list=ordered_basis2,diag_terms=kinetic_diag,tun_terms=test_tun_dic,UHK=1,mu=UHK/2)

# testf1={0:makeH.px,1:makeH.p0}
# testf2={0:makeH.px,1:makeH.p0}
# testh2=makeH.tpptwice2(state_list=mono_basis2,pauli_dic1=testf1,pauli_dic2=testf2,Uff=1,dof=2)
# testsymc=np.transpose(np.conjugate(testsym))
# prod=np.dot(testsymc,np.dot(testh2,testsym))
# print(np.allclose(prod,testh2))



####Let's just do the monolayer first
#For the F-terms

def Ftermcheck(generators_list,state_list,dof):
    invariants=[]
    paulidic={0:makeH.p0,1:makeH.px,2:makeH.py,3:makeH.pz}
    symmetrymatrices=[]
    for g in generators_list:
        symmetrymatrices.append(tpp2(state_list=state_list,pauli_dic=g[0],prefactor=g[1]))
    paulislist=list(itertools.product(paulidic.keys(),repeat=2))
    for p1 in paulislist:
        for p2 in paulislist:
            temppd1={0:paulidic[p1[0]],1:paulidic[p1[1]]}
            temppd2={0:paulidic[p2[0]],1:paulidic[p2[1]]}
            Fterm=makeH.tpptwice2(state_list=state_list,pauli_dic1=temppd1,pauli_dic2=temppd2,dof=dof,Uff=1)
            initial=True
            for g in symmetrymatrices:
                transformedF=np.dot(np.conjugate(np.transpose(g)),np.dot(Fterm,g))
                if np.allclose(transformedF,Fterm):
                    initial=True
                else:
                    initial=False
                    break
            if initial:
                invariants.append(((temppd1,temppd2),Fterm))
    return invariants
    #First need to make all the pauli dictionaries so t

def inspect_elements_mono(matrix,state_list):
    for i in range(len(state_list)):
        testwf=np.zeros(len(state_list),dtype=complex)
        testwf[i]=1
        reswf=matrix.dot(testwf)
        for j in range(len(reswf)):
            if (np.isclose(np.abs(reswf[j]),0))==False:
                #kvalue1=kdic[state_list[i].particle_dic[1].dof_dic[0]]
                #kvalue2=kdic[state_list[i].particle_dic[2].dof_dic[0]]
                #theta1=get_theta(kvalue1)-theta0/2
                #theta2=get_theta(kvalue2)-theta0/2
                #coeff1=np.cos(theta1)*v*np.sqrt(np.dot(kvalue1,kvalue1))
                #coeff2=np.cos(theta2)*v*np.sqrt(np.dot(kvalue2,kvalue2))
                print(f"Initial state: {makeH.eyepreservation_mono(state_list[i])}, Output state: {makeH.eyepreservation_mono(state_list[j])}, Amplitude {reswf[j]}")#, c1: {coeff1},c2:{coeff2}
        print(f"NEXT STATE {makeH.eyepreservation(state_list[i])}")
        #exit()


def inspect_elements(matrix,state_list):
    for i in range(len(state_list)):
        testwf=np.zeros(len(state_list),dtype=complex)
        testwf[i]=1
        reswf=matrix.dot(testwf)
        for j in range(len(reswf)):
            if (np.isclose(np.abs(reswf[j]),0))==False:
                #kvalue1=kdic[state_list[i].particle_dic[1].dof_dic[0]]
                #kvalue2=kdic[state_list[i].particle_dic[2].dof_dic[0]]
                #theta1=get_theta(kvalue1)-theta0/2
                #theta2=get_theta(kvalue2)-theta0/2
                #coeff1=np.cos(theta1)*v*np.sqrt(np.dot(kvalue1,kvalue1))
                #coeff2=np.cos(theta2)*v*np.sqrt(np.dot(kvalue2,kvalue2))
                print(f"Initial state: {makeH.eyepreservation(state_list[i])}, Output state: {makeH.eyepreservation(state_list[j])}, Amplitude {reswf[j]}")#, c1: {coeff1},c2:{coeff2}
        print(f"NEXT STATE {makeH.eyepreservation(state_list[i])}")
        #exit()


def Ntermcheck(generators_list,state_list,dof):
    invariants=[]
    paulidic={0:makeH.p0,1:makeH.px,2:makeH.py,3:makeH.pz}
    symmetrymatrices=[]
    for g in generators_list:
        symmetrymatrices.append(tpp2(state_list=state_list,pauli_dic=g[0],prefactor=g[1]))
    paulislist=list(itertools.product(paulidic.keys(),repeat=2))
    for p1 in paulislist:
        temppd1={0:paulidic[p1[0]],1:paulidic[p1[1]]}
        Nterm=makeH.HK_N_monolayer2(state_list=state_list,pd1=temppd1,Un=1)
        initial=True
        for g in symmetrymatrices:
            transformedF=np.dot(np.conjugate(np.transpose(g)),np.dot(Nterm,g))
            if np.allclose(transformedF,Nterm):
                initial=True
            else:
                initial=False
                break
        if initial:
            invariants.append((temppd1,Nterm))
    return invariants
    #First need to make all the pauli dictionaries so t

#region quick tests for N and F terms
#########################
#Test for N and F term checks
# invlist=Ntermcheck(generators_list=sym_generators,state_list=mono_basis2,dof=2)
# nzinvlist=[x for x in invlist if np.allclose(0,x[1])==False]
# #print(nzinvlist)

# # print(len(invlist))
# count=0
# for i in [x[0] for x in nzinvlist]:
#     count+=1
#     print(f"{count}. ({0}:{i[0].__name__},{1}:{i[1].__name__})")

# print(len([x[0] for x in nzinvlist]))

# testF=makeH.tpptwice2(state_list=mono_basis2,pauli_dic1={0:makeH.pz,1:makeH.pz},pauli_dic2={0:makeH.pz,1:makeH.p0},Uff=1,dof=2)
# inspect_elements_mono(state_list=mono_basis2,matrix=testF)
# print(testF)
###########################
#endregion


tqs=[np.array([-1,-1]),np.array([1,0]),np.array([0,1])]
testnonlayer={'sublattice':2,'spin':2}

basis1=makeH.generate_shell_basis(shell_count=2,q_vecs=tqs,number_of_particles=2,nonlayer=testnonlayer)

# for s in basis1:
#     print(makeH.eyepreservation(s))


testsym=tpp3(state_list=basis1,pauli_dic=c3z[0],prefactor=c3z[1])
# print(testsym)
inspect_elements(matrix=testsym,state_list=basis1)