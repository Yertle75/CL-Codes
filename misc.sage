

def dmin_crossed(c,y,c_star):
	if (c-y).hamming_weight() >= (c_star-y).hamming_weight():
		return True
	else:
		return False

def sub_vector(v,e):
	v_sub = []
	for i,b in enumerate(e):
		if b != 0:
			v_sub.append(v[i])
	return vector(v_sub)

def sub_matrix(M,e):
	M_sub = []
	cols = M.columns()
	for i,b in enumerate(e):
		if b != 0:
			M_sub.append(cols[i])
	return matrix(M_sub).transpose()


