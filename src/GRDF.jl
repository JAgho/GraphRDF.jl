module GRDF
using StaticArrays, LinearAlgebra
struct Basis{T, R} 
    vec::SVector{T, R}
end

struct Edge{T, R}
    a::SVector{T, R}
    b::SVector{T, R}
    ρ::R
end

basis(a) = (a ./ norm(a)^2)
N = 30
basis_vectors = basis.(reinterpret(SVector{2, Float64}, exp.(Complex.(0, range(0, 2*pi, N)))))
res = [zeros(100) for _ in 1:length(basis_vectors)]




function cross_conv(x::Edge, y::Edge, basis_vectors)
    for basis in basis_vectors
        x1 = x.a ⋅ basis
        x2 = x.b ⋅ basis
        y1 = y.a ⋅ basis
        y2 = y.b ⋅ basis
        trapnodes = sort([x1,x2,y1,y2])
        res .= proj_pos_a()

    end
end


# Write your package code here.

end