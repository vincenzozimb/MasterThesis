using ITensors

let 

    i = Index(2)
    j = Index(2)

    A = [1 2; 2 3]
    A = ITensor(A,i,j)
    @show sqrt(18)

    
end
