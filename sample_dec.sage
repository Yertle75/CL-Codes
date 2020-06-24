import random
import sys
from time import time
from sys import stdout
import matplotlib.pyplot as plt
load('misc.sage')

####################################
# In this decoding technique we just take a bunch of 
# random localisators and take their union, if there are some bad locators in our sample, the decoding will likely fail and we try 
# until it works

##################################
# Other functionalities

def vector_intersection(l):
	n = len(l[0])
	c = vector([1 for i in range(n)])
	for i in range(n):
		for x in l:
			c[i] = c[i] and x[i]
	return c

def vector_union(l):
	n = len(l[0])
	c = vector(GF(2),n)
	for i in range(n):
		for x in l:
			c[i] = c[i] or x[i]
	return c

def good_random_vector(k,l):
	v = vector(GF(2),k)
	l.append(v)
	while v in l:
		v = random_vector(GF(2),k)
	return v

def sample_dec(y,p,L,C):
	dim_L = L.rank()
	### Locators extraction
	S,Svs,dim_S = find_loc_space(y,L,GF(2))
	if dim_S+p == dim_L-1:
		print 'WOW THATS RARE: '+str(dim_S)
	if dim_S+p == dim_L: #easy case : all locators are good locators 
		l = vector_union(S)
		return loc_solve(y,l,C)
	else:
		for i in range(20):
			g1 = good_random_vector(dim_S,[])
			g2 = good_random_vector(dim_S,[g1])
			#g3 = good_random_vector(dim_S,[g1,g2])
			#g4 = good_random_vector(dim_S,[g1,g2,g3])
			l = vector_union([g1*S,g2*S])
			# l = vector_union([g1*S,g2*S,g3*S,g4*S])
			c_star = loc_solve(y,l,C)
			if c_star != None and (c_star+y).hamming_weight()<=p:
				return c_star
		return vector(GF(2),256)



