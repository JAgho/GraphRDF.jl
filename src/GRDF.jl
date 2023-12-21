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

# dot(point, basis) = projected point

proj(a, b::SVector{2, Float64}) = dot(a,b)

const R2 = SVector{2, Float64}

function conv_edges(a::Tuple{R2, R2}, b::Tuple{R2, R2}, r::R2; N = 100, rmax = 1, res = zeros())
    x = [a[1] ⋅ r, a[2] ⋅ r]
    y = [b[1] ⋅ r, b[2] ⋅ r]
    x_left = x[1] < x[2] ? (1, 2) : (2, 1) # get sort indices for dot products
    y_left = y[1] < y[2] ? (1, 2) : (2, 1) 
    sx = x[x_left] # create temporary sorted vars
    sy = y[y_left]
    sl = sx[1] - sx[2] > sy[1] - sy[2] ? (sy, sx) : (sx, sy) #order from smallest top hat to largest
    return sl[2][1] - sl[1][2], sl[2][1] - sl[1][2], sl[2][2] - sl[1][1]
end

function trap2hist(trap, h, rmax, N, res=zeros(N))
    grad = h/(trap[1] - trap[2])
    l1 = floor(trap[1])
end




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