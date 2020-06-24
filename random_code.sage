for i in range(16):
     L = matrix(QAryReedMullerCode(GF(16),i,2).basis())
     LL = matrix(GF(16),[u.pairwise_product(v) for u in L for v in L])
     LL = LL.echelon_form()[:LL.rank()]
     LL_ = binary_dual(16,LL)
     print(L.rank(),LL.rank(),256-LL_.rank())

L = matrix(QAryReedMullerCode(GF(16),i,2).basis())
LL = matrix(GF(16),[u.pairwise_product(v) for u in L for v in L])
LL = LL.echelon_form()[:LL.rank()]
LL_ = binary_dual(16,LL)
V16 = VectorSpace(GF(16),256)
V2 = VectorSpace(GF(2),256)
LLvs = V16.subspace_with_basis(LL)
LL_vs = V2.subspace_with_basis(LL_)
C = matrix(LLvs.complement().basis())
C_ = matrix(LL_vs.complement().basis())

print(L.rank(),LL.rank(),C.rank(),C_.rank())

S,Svs,dim_S = find_loc_space(y,C)
cdict = gen_dict(S)
plt.plot(cdict.values()) ; plt.show()












def find_error_vector(p,S):
	cols = S.columns()
	small = cols[:p]
	M = matrix(GF(2),small)
	r = rank(M)
	for i in range(p,len(cols)):
		for j in range(p):
			save = subset[j]
			subset[j] = cols[i]
			M = matrix(GF(2),subset)
			if M.rank()<








### print loc repartition 
import matplotlib.pyplot as plt

values = [0 for i in range(201)]

for i in range(N):
	c,y,e = random_codeword(p,C)
	S,Svs = find_loc_space(y,C)
	for l in Svs[1:]: 
		hw = l.hamming_weight()
		values[hw]+=1

plt.plot(values)
plt.show()

def random_code_n_choose_k(n,t,k1,k2):
	print('Generating code...')
	### create the set of t choose k possible columns
	positions = list(itertools.combinations(range(t), k1))
	positions.extend(list(itertools.combinations(range(t), k2)))
	kept_positions = random.sample(positions,n)
	L = matrix(GF(2),t,0)
	for pos in kept_positions:
		v = vector(GF(2),t)
		for i in pos:
			v[i]=1
		L = L.augment(v)
	print('| Generated L : dim = '+str(L.rank()))
	L.echelonize()
	return L[:L.rank()]



def skim_loc_space(y,S,C):
	S_skimmed = []
	for l in S:
		flag,c = loc_solve(y,l,C)
		if flag:
			S_skimmed.append(l)
	return matrix(GF(2),S_skimmed)

def find_majority(k):
    myMap = {}
    maximum = ( '', 0 ) # (occurring element, occurrences)
    for n in k:
        if str(n) in myMap: myMap[str(n)] += 1
        else: myMap[str(n)] = 1

        # Keep track of maximum on the go
        if myMap[str(n)] > maximum[1]: maximum = (n,myMap[str(n)])
	return maximum


def test(p,C,L):
	dim = 0
	while dim != 6:
		c,y,e = random_codeword(p,C)
		S,Svs = find_loc_space(y,C)
		dim = Svs.dimension()
		print dim
	_,_,e_fake = random_codeword(p,C)
	red_L = []
	cols = L.columns()
	for i,b in enumerate(e):
		if b== 1:
			red_L.append(cols[i])
	red_L = matrix(GF(2),red_L)
	print(red_L.transpose().echelon_form())
	bad = []
	for l in Svs:
		bad.append(l)
	#inter_bad = vector_intersection(bad)
	cols = matrix(GF(2),bad).columns()
	red_inter_bad = []
	for i,b in enumerate(e):
		if b== 1:
			red_inter_bad.append(cols[i])
	red_inter_bad = matrix(GF(2),red_inter_bad)
	print 'bad'
	print red_inter_bad.transpose()
	print 'y'
	red_L = []
	cols = y
	for i,b in enumerate(e):
		if b== 1:
			red_L.append(cols[i])
	print(matrix(GF(2),red_L))
	return red_L


def decode_v2(c,y,e,L,C,logs):
	n = len(C[0])
	zeros = vector(GF(2),n)
	### Locators extraction
	S,Svs = find_loc_space(y,C)
	S_ = S
	### let's first skim our space of obvious bad locators
	print 'loc dim before skimming: '+str(S.rank())
	S = skim_loc_space(y,S,C)
	print 'loc dim after : '+str(S.rank())
	l = vector_union(S)
	lC = matrix(GF(2),[l.pairwise_product(u) for u in C])
	yl = y.pairwise_product(l)
	try:
		g = lC.solve_left(yl)
	except:
		print('solve error')
		zeros = vector(GF(2),len(C[0]))
		print ' in S_skimmed'
		for l in Svs:
			inter = vector_intersection([e,l])
			print inter == zeros
		print ' in S'
		return 0
		
	K = lC.kernel()
	if K.rank() > 0:
		print 'kernel error, rank: '+str(K.rank())
		return 0
	if g*C == c:
		return 1
	else:
		print 'FAAAAAAAAAAAAAAAIL'
		zeros = vector(GF(2),len(C[0]))
		for l in list(S):
			inter = vector_intersection([e,l])
			print inter == zeros
		return 0
