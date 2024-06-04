using OMEinsum, Test, LinearAlgebra, Yao

struct OMPS{T}
    tensors::Vector{Array{T, 3}}
end

function truncated_svd(tmat::AbstractMatrix, dmax::Int, atol::Real)
    u, s, v = LinearAlgebra.svd(tmat)
    dmax = min(searchsortedfirst(s, atol, rev=true), dmax, length(s))
    return u[:, 1:dmax], s[1:dmax], v'[1:dmax, :], sum(s[dmax+1:end])
end

function center_canonicalize!(mps::OMPS, l::Int) # l is the orthogonality center
    for i = 1:l-1 # canonicalize left part
        left, right = mps.tensors[i], mps.tensors[i+1]
        mleft, mright = reshape(left, :, size(left, 3)), reshape(right, size(right, 1), :)
        u, s, v, trunc = truncated_svd(mleft * mright, size(mleft, 2), 0)
        mps.tensors[i] = reshape(u, size(left, 1), size(left, 2), size(u, 2))
        mps.tensors[i+1] = reshape(Diagonal(s) * v, size(v, 1), size(right, 2), size(right, 3))
        mps.tensors[i+1] = permutedims(mps.tensors[i+1], (2, 1, 3))
    end

    for i = length(mps.tensors)-1:-1:l # canonicalize right part
        left, right = mps.tensors[i], mps.tensors[i+1]
        mleft, mright = reshape(left, :, size(left, 3)), reshape(right, size(right, 1), :)
        u, s, v, trunc = truncated_svd(mleft * mright, size(mleft, 2), 0)
        mps.tensors[i] = reshape(u * Diagonal(s), size(left, 1), size(left, 2), size(u, 2))
        mps.tensors[i+1] = reshape(v, size(v, 1), size(right, 2), size(right, 3))
        mps.tensors[i+1] = permutedims(mps.tensors[i+1], (2, 1, 3))
    end
    return mps
end
#spin-up is |0>, spin-down is |1>



t1 = rand(ComplexF64, 1, 2, 2)
t2 = rand(ComplexF64, 2, 2, 2)
t3 = rand(ComplexF64, 2, 2, 1)
C = [t1, t2, t3]
CC = OMPS(C)

psi = [1;
       0]

p = ein"aij,jpk,klb,ain,nqm,mlb,p,q->"(t1,t2,t3,conj.(t1),conj.(t2),conj.(t3),psi,conj.(psi)) ./ ein"aij,jpk,klb,ain,npm,mlb->"(t1,t2,t3,conj.(t1),conj.(t2),conj.(t3))

center_canonicalize!(CC,2)

#for i = 1:length(CC.tensors)
#    @show CC.tensors[i]
#end



p_can = ein"i,j,ikl,jkl->"(psi,conj.(psi),CC.tensors[2],conj.(CC.tensors[2])) ./ ein"ijk,ijk->"(CC.tensors[2],conj.(CC.tensors[2]))

#=
println( ein"ijk,ijl->kl"(CC.tensors[1],conj.(CC.tensors[1])) )
println( ein"kij,lij->kl"(CC.tensors[3],conj.(CC.tensors[3])) )
=#
@show p
@show p_can