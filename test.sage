#########
# What I usually do is generate L,C,PI in a sage console and then run load('test.sage')

from time import time
#####################################
### Initial parameters and loading

load('gen.sage')
load('lowrank_dec.sage')
load('sample_dec.sage')

##############################
### Additional parameters

p = 50
N = 1000

####################
### Test Phase

start_time = int(time())
nb_success = 0
nb_okfailure = 0
for i in range(N):
	c,y,e = random_codeword(p,C,GF(2))
	c_star = lowrank_dec(y,p,L,C,GF(2))
	if c_star == c:
		nb_success+= 1
	elif c_star == None:
		continue
	else:
		if (c_star-y).hamming_weight()<=(c-y).hamming_weight(): 
			# this counts the number of times the closest codeword is not the initial one
			nb_okfailure += 1
		
finish_time = int(time())
print('=======================')
print(' Final Result on '+str(N)+' tests in '+str(finish_time-start_time)+'s with weight '+ str(p))
print('Success = '+str(nb_success)+', OK failure = '+str(nb_okfailure))

