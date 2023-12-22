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
N = 9
basis_vectors = basis.(reinterpret(SVector{2, Float64}, exp.(Complex.(0, range(0, pi, N)))))
res = [zeros(100) for _ in 1:length(basis_vectors)]

# dot(point, basis) = projected point

proj(a, b::SVector{2, Float64}) = dot(a,b)

const R2 = SVector{2, Float64}

atest = R2((0,0))
btest = R2((0,1))
ctest = R2((0,3))
dtest = R2((0,4))

e1 = (atest, btest)
e2 = (ctest, dtest)
bas = basis_vectors[5]


fig, ax = plot([e1...])
plot!(ax, [e2...])
plot!(ax, basis_vectors[5])

@code_warntype conv_edges(e1, e2, bas)
@benchmark conv_edges(e1, e2, bas)

function sort2!(x)
    swap(x) = x[2], x[1]
    x[1]>x[2] ? swap(x) : x
end

function conv_edges(a::Tuple{R2, R2}, b::Tuple{R2, R2}, r::R2)
    x = (a[1] ⋅ r, a[2] ⋅ r)
    y = (b[1] ⋅ r, b[2] ⋅ r)
    #print("x vectors are $x, y vectors are $y")
    sx1, sx2 = sort2!(x)
    sy1, sy2 = sort2!(y)
    #order from smallest top hat to largest
    if sx1 - sx2 > sy1 - sy2
        return sy1 - sx2, sy1 - sx1, sy2 - sx2, sy2 - sx1
    else
        return sx1 - sy2, sx1 - sy1, sx2 - sy2, sx2 - sy1
    end
end

function trap2hist(trap, h, rmax, N, res=zeros(N))
    grad = h/(trap[1] - trap[2])
    ls = floor.(trap)
    rs = ls .+ 1
    leftmost = 1

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