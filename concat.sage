#####################
# Testing the ideas with concatenated code
# (only generates L)

load('all.sage')
A = codes.GeneralizedReedSolomonCode(GF(64).list(),10).basis()
B = [[0 for i in range(64*32)] for j in range(6*10)]
M = reedmuller(1,5,GF(2))
z6 = GF(64).primitive_element()
for i in range(10):
	for l in range(6):
		m = z6**l
		for j in range(64):
			a = A[i][j]*m
			b = a.integer_representation()
			b = Integer(b).digits(2)
			v = vector(GF(2),6).list()
			for k,bit in enumerate(b):
				v[k]= bit
			b = vector(GF(2),v)
			#beware b's bits are 'reversed'... but we dont care
			b = b*M
			for k in range(32):
				B[6*i+l][32*j+k] = b[k]
L = matrix(GF(2),B)		
