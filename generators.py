import numpy as np
# import makeH
# import pathdiag
# import numpy as np
# import sys


# #Global variables
# from variables import *








def generate_matrices(pairs):
    # Convert matrices to tuples for comparison
    def matrix_to_tuple(matrix):
        return tuple(map(tuple, matrix))

    # Initialize dictionary to keep track of unique matrices
    unique_matrices = {}

    # Initial population of the dictionary with the given pairs
    for string, matrix in pairs:
        mat_tuple = matrix_to_tuple(matrix)
        if mat_tuple not in unique_matrices:
            unique_matrices[mat_tuple] = string

    # Iterate over all combinations
    length = len(pairs)
    for i in range(length):
        for j in range(i, length):
            string1, matrix1 = pairs[i]
            string2, matrix2 = pairs[j]
            new_matrix = np.dot(matrix1, matrix2)
            new_string = string1 + string2

            mat_tuple = matrix_to_tuple(new_matrix)
            if mat_tuple not in unique_matrices:
                unique_matrices[mat_tuple] = new_string

    return unique_matrices


# Example usage
pairs = [
    ("X", np.array([[0, 1], [1, 0]])),  # Pauli X
    ("Y", np.array([[0, -1j], [1j, 0]])), # Pauli Y
    # Add more pairs as needed
]

# result = generate_matrices(pairs)
# for matrix, string in result.items():
#     print(f"{string}: {np.array(matrix)}")

def generate_triples(op_name,pauli_list,coeff,conjugation):
    tp_matrix=pauli_list[0]
    
    for p in range(1,len(pauli_list)):
        tp_matrix=np.kron(tp_matrix,pauli_list[p])
    tp_matrix=coeff*tp_matrix
    def pauli_func(target_matrix,target_conjugation):
        if conjugation:
            resulting_matrix=np.dot(tp_matrix,np.conjugate(target_matrix))
            resulting_conj=(conjugation!=target_conjugation)
        else:
            resulting_matrix=np.dot(tp_matrix,target_matrix)
            resulting_conj=(conjugation!=target_conjugation)
        return resulting_matrix,resulting_conj
        
    return (op_name,pauli_func,tp_matrix,conjugation)

s0=np.array([[1,0],[0,1]],dtype=complex)
sx=np.array([[0,1],[1,0]],dtype=complex)
sy=np.array([[0,-1j],[1j,0]],dtype=complex)
sz=np.array([[1,0],[0,-1]],dtype=complex)



# test=generate_triples(op_name='C3z',pauli_list=[sz,sz],coeff=1,conjugation=True)
# print(test)
#print(test[1](target_matrix=np.identity(n=2**2,dtype=complex),target_conjugation=True))



def my_generate(generator_triples):
    
    gen_triples=[]
    #1. generate pairs
    for string,paulis,coeff,conj in generator_triples:
        gen_triples.append(generate_triples(op_name=string,pauli_list=paulis,coeff=coeff,conjugation=conj))
    
    #2. Starting with the identity, multiply everything in unique matrices before and after by the generators and test whether the ending set is the same as the starting set
    dim=2**len(paulis)
    unique_matrices={'1':(np.identity(n=dim,dtype=np.complex128),False)}
    converged=False
    count=0
    while converged==False:
        print(f"count: {count}")
        current_pairs=list(unique_matrices.values())
        new_pairs=[]
        for operator_name in unique_matrices.keys():
            op_matrix=unique_matrices[operator_name][0]
            op_conj=unique_matrices[operator_name][1]
            for op_name,matrix_func,tp_matrix,conjugation in gen_triples:
                premult_pair=matrix_func(target_matrix=op_matrix,target_conjugation=op_conj)
                if check_in(premult_pair,current_pairs)==False:
                    new_pairs.append((op_name,operator_name,premult_pair))
                if op_conj:
                    postmult_pair=(np.dot(op_matrix,np.conjugate(tp_matrix)),op_conj!=conjugation)
                else:
                    postmult_pair=(np.dot(op_matrix,tp_matrix),op_conj!=conjugation)
                if check_in(postmult_pair,[premult_pair])==False and check_in(postmult_pair,current_pairs)==False:
                    new_pairs.append((operator_name,op_name,postmult_pair))
        print(new_pairs)
        if len(new_pairs)==0:
            converged=True
        else:
            for t in new_pairs:
                unique_matrices[t[0]+t[1]]=t[2]
    return unique_matrices
                


        

def check_in(pair,list):
    pmat=pair[0]
    pconj=pair[1]
    inlist=False
    for list_pair in list:
        mat_match=np.allclose(np.real(pmat),np.real(list_pair[0])) and np.allclose(np.imag(pmat),np.imag(list_pair[0]))
        #print(f'mat_match {mat_match}')
        conj_match=(pconj==list_pair[1])
        #print(f'conj_match {conj_match}')
        #print(pconj)
        
        if mat_match and conj_match:
            inlist=True
            break
    return inlist


# test_pair=(np.kron(),False)
# list1=[(sy,True),(sx,False)]
# print(check_in(test_pair,list1))

#test_triples_list=[generate_triples(op_name='C3z',pauli_list=[sz,sz],coeff=1,conjugation=True)]
test_triples_list=[['C3z',[s0,sz],1j,False],['C2z',[sx,s0],1,False]]

print(test_triples_list)
test=my_generate(generator_triples=test_triples_list)
print(test.keys())

    


    



