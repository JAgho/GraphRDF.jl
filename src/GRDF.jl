using StaticArrays, LinearAlgebra

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

#discretise a linear function with compact support onto results vector res.
#the function is bounded lx < x < rx and spans from ly->ry
#this definitely needs functionalisation and optimisation as binwidth and res are global
function piecewise_int!(lx, rx, ly, ry, res)
    g = ry-ly / abs(rx - lx)
    m = abs(g)
    N = length(res)
    ri = bin_offset(rx)
    li = bin_offset(lx)
    re = bin_index(ri)
    le = bin_index(li)
    if le+1 > N
        return
    end
    rloc = (ri) * binwidth
    lloc = (li+1) * binwidth
    if le+1 < re
        w1 = abs(lx -lloc)
        lyl = min(ly + g*w1, ly)
        area1 = w1*lyl + (w1*w1*m/2)
        w2 = abs(rx - rloc)
        ryl = min(ry - g*w2, ry)
        area2 = w2*ryl + (w2*w2*m/2)
        res[le] += area1
        res[re] += area2
        area_increment =  g*binwidth*binwidth
        area = ly*binwidth + area_increment
        for i in le+1:min(re-1, N)
            res[i] += area
            area += area_increment
        end
    elseif le < re
        w1 = abs(lx -lloc)
        lyl = min(ly + g*w1, ly)
        area1 = w1*lyl + (w1*w1*m/2)
        w2 = abs(rx - rloc)
        ryl = min(ry - g*w2, ry)
        area2 = w2*ryl + (w2*w2*m/2)
        res[le] += area1
        res[re] += area2
    else
        res[le] += m*abs(rx - lx)/2
    end
end

#computes the discretisation of a trapezoid onto results vector res.
#trapezoid is a 4-tuple and height is z
function piecewise_trap!(trap, z, res)
    piecewise_int!(trap[1], trap[2], 0, z, res)
    piecewise_int!(trap[2], trap[3], z, z, res)
    piecewise_int!(trap[3], trap[4], z, 0, res)
end


N = 9
basis_vectors = basis.(reinterpret(SVector{2, Float64}, exp.(Complex.(0, range(0, pi, N)))))
rdf = [zeros(100) for _ in 1:length(basis_vectors)]

atest = R2((0,0))
btest = R2((0,1))
ctest = R2((0,3))
dtest = R2((0,5))

e1 = (atest, btest)
e2 = (ctest, dtest)
bas = basis_vectors[5]

trap = conv_edges(e2, e1, bas)
nbins = 101 #no of bins we work over
res=zeros(nbins) #preallocated array for results
rmax = 6.0 #maximal distance

mid = 50
binwidth = rmax/floor(nbins)*2
pos = 1.0
bin_offset(pos) = floor(Int64, pos / binwidth)
bin_index(offset) = offset + mid


piecewise_trap!(trap, z, res)


