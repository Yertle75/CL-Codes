import itertools
import random
from sage.coding.reed_muller_code import QAryReedMullerCode
from sage.coding.reed_muller_code import BinaryReedMullerCode

#########################
# Other functionalities

def random_codeword(p,C,K):
	n = len(C[0])
	k = C.rank()   # code dimension
	g = random_vector(K,k)  # generator
	c = g * C   #codeword
	e = vector(K,n)  # error vector
	for i in sample(xrange(n),p):
		e[i]=1
	y=c+e    #errored codeword
	return c,y,e




######################
# random codes (experimented with before Reed-Muller )

def add_ones(L,K):
	n = len(L[0])
	#add a row of ones to the matrix
	ones = vector(K,n,[1 for i in range(n)])
	L = L.transpose().augment(ones).transpose()   # uh !
	print('| Generated L	: '+str(rank(L))+' x '+str(n))
	return L.echelon_form()[:L.rank()]

def random_code(nb_lines,nb_cols,K = _K_):
	print('Generating code ...')
	V = VectorSpace(K,nb_lines)
	R = random.sample(list(V)[1:],nb_cols)
	L = matrix(K,nb_lines,0)
	for v in R:
		L = L.augment(v)
	print('| Generated L	: '+str(nb_lines)+' x '+str(nb_cols))
	return L.echelon_form()[:L.rank()]

#########################
# Reed-Muller codes 

def reedmuller(order,nb_vars,K):
	print('Generating code ...')
	if K == GF(2):
		L = matrix(BinaryReedMullerCode(order,nb_vars).basis())
	else:
		L = matrix(QAryReedMullerCode(K, order, nb_vars).basis())
	n = len(L[0])
	k = rank(L)
	print('| Generated L	: '+str(k)+' x '+str(n))
	return L

#################################		
# generate the LL, C and PI spaces

def generate_code_family(L,K):
	print('Generating code family...')
	n = len(L[0])
	V = VectorSpace(K,n)
	LL = matrix(K,[u.pairwise_product(v) for u in L for v in L])
	LL = LL.echelon_form()[:LL.rank()]
	print('| Generated LL	: '+str(LL.rank())+' x '+str(n))
	LLvs = V.subspace_with_basis(LL)
	C = matrix(LLvs.complement().basis())
	print('| Generated C	: '+str(C.rank())+' x '+str(n))
	PI = matrix(K,[u.pairwise_product(v) for u in L for v in C])
	PI = PI.echelon_form()[:PI.rank()]
	print('| Generated PI	: '+str(PI.rank())+' x '+str(n))
	return LL,C,PI
##################################
# generate the LL, and binary C spaces 
#
# This one is used specifically when L is qary and we want to extract the binary subcode of C
#
def binary_equivalent(M,K):
	n = len(M[0])
	B = matrix(GF(2),n,0)
	for l in M:
		p = [x.polynomial() for x in l]
		for i in range(K.order()):
			b = [x[i] for x in p]
			B = B.augment(vector(GF(2),b))
	return B.transpose().echelon_form()[:B.rank()]


def generate_binary_family(L,K):
	print('Generating binary code family...')
	n = len(L[0])
	Vbin = VectorSpace(GF(2),n)
	LL = matrix(K,[u.pairwise_product(v) for u in L for v in L])
	LL = LL.echelon_form()[:LL.rank()]
	LLbin = binary_equivalent(LL,K)
	print('| Generated LL	: '+str(LL.rank())+' x '+str(n))
	LLbinvs = Vbin.subspace_with_basis(LLbin)
	C = matrix(LLbinvs.complement().basis())
	print('| Generated C	: '+str(C.rank())+' x '+str(n))
	PI = matrix(K,[u.pairwise_product(v) for u in L for v in C])
	PI = PI.echelon_form()[:PI.rank()]
	print('| Generated PI	: '+str(PI.rank())+' x '+str(n))
	return LL,C,PI
