using ITensors

println("\n\n")

i = Index(2)
j = Index(5)
k = Index(2)

A = randomITensor(i, j, k)
@info inds(A)

U, S, V = svd(A, (i, k))

@info inds(U)
# @info inds(S)
# @info inds(V)
