##########################
# The basic functionalities are: 
# find_loc_space : given y find S the set of locators for that y
# loc_solve : given y,l find a c_star such that c_star * l = y * l
# exhaustive_kernel_search : given y and a set of potential c_star, explore them all to find the closest to y

#####################
# lowrank_dec uses the following heuristic : knowing the rank of S|e, we know the number of columns we need to isolate.
# We then sort all the columns by appearance ( beware that we are talking of the vectorial space spanned by columns so if two columns
# are taken, their sum must be too) and try to solve that way.

def find_loc_space(y,L, K):
	### Locators extraction
	n = len(L[0])
	yL = matrix([u.pairwise_product(y) for u in L])
	R = yL*L.transpose()
	R_ker = matrix(R.kernel().basis())
	
	S = R_ker*L  # S is the set of all possible locators
	Svs = VectorSpace(K,n).subspace_with_basis(S)
	dim_S = Svs.dimension()
	return S,Svs,dim_S

def cweight(a):
	(_,x) = a
	return x

def cweight2(a):
	(_,x),(_,y) = a
	return (max(x,y),min(x,y))

def cweight3(a):
	(_,x),(_,y),(_,z) = a
	l = [x,y,z]
	l.sort(reverse = True)
	return (l) #pas tout a fait precis mais ok


def gen_dict(S):
	cdict = {}
	for c in S.columns():
		c = tuple(c)
		if c in cdict: cdict[c] += 1
		else: cdict[c] = 1
	return cdict

def find_max(cdict):
	k_star = ''
	v_star = 0
	for k,v in cdict.items():
		if v>v_star:
			k_star,v_star = k,v
	cdict.pop(k_star)
	return k_star,v_star		

def improve_loc(l,c,S):
	for i,c_ in enumerate(S.columns()):
		if c == c_:
			l[i] = 0


### Test on random error vector
def exhaustive_kernel_search(g,y,R_kernel,C,K):
	n = len(C[0])
	c_star = vector(K,n)
	for v in R_kernel:
		c_guess = (g+v)*C
		if (c_guess+y).hamming_weight() < (c_star+y).hamming_weight():
			c_star = c_guess
	return c_star
		
def loc_solve(y,l,C,K):
	lC = sub_matrix(C,l)
	yl = sub_vector(y,l)
	try:
		g = lC.solve_left(yl)
	except:
		return None # This locator (or loc_union) was bad
	c_star = g*C
	lC_kernel = lC.kernel()
	ker_dim = lC_kernel.dimension()
	if ker_dim > 5:
		### kernel too big for exhaustive search, fingers crossed, we pick one at random
		print 'warning big kernel : dim K = '+str(ker_dim)
		return c_star
	elif ker_dim > 0:
		### we hope to get better results by looking the closest vect in the kernel
		#print 'warning intermediate kernel : dim K = '+str(ker_dim)
		c_star = exhaustive_kernel_search(g,y,lC_kernel,C,K)
		return c_star
	else:
		### only one solution, good
		return c_star


def lowrank_dec(y,p,L,C,K):
	n = len(L[0])
	dim_L = rank(L)
	### Locators extraction
	S,Svs,dim_S = find_loc_space(y,L,K)
	### generate a dict containing columns
	cdict = gen_dict(S)
	
	###generate the basic locator removing zeros
	l = vector(K,[1 for i in range(n)])
	zeros = vector(K,dim_S)
	for i,c in enumerate(S.columns()):
		if c == zeros:
			l[i] = 0
			
	### different cases
	if dim_S+p == dim_L:
		c_star = loc_solve(y,l,C,K)
		return c_star
	elif dim_S+p == dim_L+2:
		if tuple(zeros) in cdict:
			cdict.pop(tuple(zeros))
		clist = cdict.items()
		clist.sort(key = cweight, reverse = True)
		for (c,v) in clist:
			l_ = copy(l)
			improve_loc(l_,vector(K,c),S)
			c_star = loc_solve(y,l_,C,K)
			if c_star != None and (c_star-y).hamming_weight() <= p:
				return c_star
	elif dim_S+p == dim_L+4:
		if tuple(zeros) in cdict:
			cdict.pop(tuple(zeros))
		clist = cdict.items()
		clist = list(itertools.combinations(clist, 2))
		clist.sort(key = cweight2, reverse = True)
		for (c1,v1),(c2,v2) in clist:
			l_ = copy(l)
			improve_loc(l_,vector(K,c1),S)
			improve_loc(l_,vector(K,c2),S)
			improve_loc(l_,vector(K,c1)+vector(K,c2),S)
			c_star = loc_solve(y,l_,C,K)
			if c_star != None and (c_star-y).hamming_weight() <= p:
				return c_star
	elif dim_S+p == dim_L+6:
		if tuple(zeros) in cdict:
			cdict.pop(tuple(zeros))
		clist = cdict.items()
		clist = list(itertools.combinations(clist, 3))
		clist.sort(key = cweight3, reverse = True)
		for (c1,v1),(c2,v2),(c3,v3) in clist:
			l_ = copy(l)
			improve_loc(l_,vector(K,c1),S)
			improve_loc(l_,vector(K,c2),S)
			improve_loc(l_,vector(K,c3),S)
			improve_loc(l_,vector(K,c1)+vector(K,c2),S)
			improve_loc(l_,vector(K,c1)+vector(K,c3),S)
			improve_loc(l_,vector(K,c2)+vector(K,c3),S)
			improve_loc(l_,vector(K,c1)+vector(K,c2)+vector(K,c3),S)
			c_star = loc_solve(y,l_,C,K)
			if c_star != None and (c_star-y).hamming_weight() <= p:
				return c_star
	else:
		print 'unsupported case below'
	return vector(K,_n_)
	
