#############
# generate a code L from projective geometry design

PG = designs.ProjectiveGeometryDesign(3, 2, GF(8))
L = matrix(GF(2),PG.incidence_matrix())
L = L.echelon_form()[:rank(L)]


