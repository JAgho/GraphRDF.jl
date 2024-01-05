using StaticArrays, LinearAlgebra

const R2 = SVector{2, Float64}
const R3 = SVector{3, Float64}

struct Edge
    e::Tuple{R2, R2}
    rho::Float64
end

basis(a) = (a ./ norm(a)^2)
proj(a, b::SVector{2, Float64}) = dot(a,b)

const R2 = SVector{2, Float64}
const R3 = SVector{3, Float64}

#fast sorting of a 2-parameter tuple / vector
function sort2!(x)
    swap(x) = x[2], x[1]
    x[1]>x[2] ? swap(x) : x
end

#convolves a pair of top hat functions with densities ρ1 and ρ2, returning  
# the y height of the trapezoid and the x nodes of its inflections
function conv_edges(a::Edge, b::Edge, r::R2)
    ae1 = a.e[1]
    ae2 = a.e[2]
    be1 = b.e[1]
    be2 = b.e[2]
    x = (ae1 ⋅ r, ae2 ⋅ r)
    y = (be1 ⋅ r, be2 ⋅ r)
    ρ = (norm(ae1 - ae2)/ norm(x[1] - x[2])) * (norm(be1 - be2) / norm(y[1] - y[2])) * a.rho * b.rho
    #print("x vectors are $x, y vectors are $y")
    sx1, sx2 = sort2!(x)
    sy1, sy2 = sort2!(y)
    #order from smallest top hat to largest
    if sx1 - sx2 > sy1 - sy2
        return (sy1 - sx2, sy1 - sx1, sy2 - sx2, sy2 - sx1), ρ
    else
        return (sx1 - sy2, sx1 - sy1, sx2 - sy2, sx2 - sy1), ρ
    end
end

#discretise a linear function with compact support onto results vector res.
#the function is bounded lx < x < rx and spans from ly->ry
#this definitely needs functionalisation and optimisation as binwidth and res are global
function piecewise_int!(lx, rx, ly, ry, res, binwidth)
    g = (ry-ly) / abs(rx - lx)
    # @show g
    # @show ry-ly
    # @show rx-lx
    m = abs(g)
    N = length(res)
    ri = bin_offset(rx)
    li = bin_offset(lx)
    re = bin_index(ri)
    le = bin_index(li)
    #@show le, re
    if le+1 > N
        return nothing
    end
    rloc = (ri) * binwidth
    lloc = (li+1) * binwidth
    if le+1 < re
        w1 = abs(lx -lloc)
        lyl = min(ly + g*w1, ly)
        area1 = w1*lyl + (w1*w1*m*0.5)
        w2 = abs(rx - rloc)
        ryl = min(ry - g*w2, ry)
        area2 = w2*ryl + (w2*w2*m*0.5)
        @inbounds res[le] += area1
        if re < N
            @inbounds res[re] += area2
        end
        area_increment =  g*binwidth*binwidth
        area = ly*binwidth + area_increment
        @inbounds for i in le+1:min(re-1, N)
            res[i] += area
            area += area_increment
        end
        @show "bwa"
    elseif le < re
        w1 = abs(lx -lloc)
        lyl = min(ly + g*w1, ly)
        area1 = w1*lyl + (w1*w1*m*0.5)
        w2 = abs(rx - rloc)
        ryl = min(ry - g*w2, ry)
        area2 = w2*ryl + (w2*w2*m*0.5)
        @inbounds res[le] += area1
        @inbounds res[re] += area2
        @show "hello"
    else
        @show "hi"
        @inbounds res[le] += m*abs(rx - lx)*0.5
    end
    return nothing
end

#computes the discretisation of a trapezoid onto results vector res.
#trapezoid is a 4-tuple and height is z
function piecewise_trap!(trap, z, res, binwidth)
    piecewise_int!(trap[1], trap[2], 0, z, res, binwidth)
    piecewise_int!(trap[2], trap[3], z, z, res, binwidth)
    piecewise_int!(trap[3], trap[4], z, 0, res, binwidth)
end

function radial_cross_corr(e1::Edge, e2::Edge, basis_vectors, rdf)
    for (result, basis) in  zip(rdf, basis_vectors)
        trap1, ρ = conv_edges(e2, e1, basis)
        piecewise_trap!(trap1, ρ, result, binwidth)
    end
end

function plot_rdf(nt, rdf)
    nt = 0
    begin
        fig, ax = lines(rdf[nt])
        for i in nt:N-nt
            lines!(ax, rdf[i])
        end
        display(fig)
    end
end

function basis_2d(N, nbins)
    basis_vectors = basis.(reinterpret(SVector{2, Float64}, exp.(Complex.(0, range(0, pi/2, N)))))
    rdf = [zeros(nbins) for _ in 1:N]
    return basis_vectors, rdf
end

bin_offset(pos)::Int64 = floor(Int64, pos / binwidth)
bin_index(offset)::Int64 = offset + mid

N = 9
nbins = 101
# basis_vectors = basis.(reinterpret(SVector{2, Float64}, exp.(Complex.(0, range(0, pi/2, N)))))
# rdf = [zeros(101) for _ in 1:length(basis_vectors)]
basis_vectors, rdf = basis_2d(N, nbins)

rmax = 6.0 #maximal distance
mid::Int64 = 50
binwidth::Float64 = rmax/floor(nbins)*2

atest = R2((0,0))
btest = R2((0,1))
ctest = R2((0,3))
dtest = R2((0,5))

e1 = (atest, btest)
e2 = (ctest, dtest)
bas = basis_vectors[5]

edge1 = Edge(e1, 1)
edge2 = Edge(e2, 1)

radial_cross_corr(edge1, edge2, basis_vectors, rdf)



# res .= 0
# ttrap = (2.0,3.0,4.0,5.0)
# ttrap = (1.0,2.0,3.0,4.0)
# ttrap = (1.8477590650225733, 2.77163859753386, 3.695518130045147, 4.619397662556434)
# ttrap = (0.13080625846028615, 0.19620938769042923, 0.2616125169205723, 0.3270156461507154)
# ttrap = (0.025182441996913264, 0.0377736629953699, 0.050364883993826534, 0.06295610499228317)
# ttrap = (0.03777241530507665, 0.05665862295761497, 0.0755448306101533, 0.09443103826269161)
# piecewise_trap!(ttrap, 1.0, res, binwidth)
# lines(res)
#plot(res)