using Oscar
using Plots

mm(A, w) = x -> sum([w[i]*abs(x[i])*A[:,i] for i in eachindex(w)])
param(A) = x -> [prod([x[i]^A[i, j] for i in eachindex(x)]) for j in eachindex(A[1,:])]
proj_mm(A, w) = x -> sum([w[i]*abs(x[i])*A[:,i] for i in eachindex(w)])/(sum([abs(w[i]*x[i]) for i in eachindex(x)]))

function sample_square(α, β, d)
    return [rand(α:.0001:β, 1, d) for _ in 1:10000]
end

function sample_square(α, d)
    return sample_square(0, α, d)
end

sample_square(3, 2)

function map_square(plt, α, d::Int, f)
    plt_pts = hcat([f(p) for p in sample_square(α, d)]...)
    scatter!(plt, plt_pts[1,:], plt_pts[2,:], label="$(α)")
end

function map_square(plt, α, A, w)
    d, n = size(A)
    ψ = mm(A, w)
    ϕ = param(A)
    map_square(plt, α, d, ψ ∘ ϕ)
end

function proj_map_square(plt, α, A, w)
    d, n = size(A)
    ψ = proj_mm(A, w)
    ϕ = param(A)
    sample_pts = vcat([vcat(1, p...) for p in sample_square()])
end

##############
## Example
##############

## from the book:
m1(x) = [x[1]*x[2] + x[1]*x[2]^2 + 2*x[2]*x[1]^2, x[1]*x[2] + 2*x[1]*x[2]^2 + x[2]*x[1]^2]

A = hcat([[1, 1], [1, 2], [2, 1]]...)
Y = toric_ideal(transpose(A))

## parametrization + moment map
plt = plot(xlims=[0,4], ylims=[0,4])
wt = [1,1,1]
map_square(plt, 3, A, wt)
map_square(plt, 2, A, wt)
map_square(plt, 1, A, wt)

wt = [100,1,1]
map_square(plt, 3, A, wt)
map_square(plt, 2, A, wt)
map_square(plt, 1, A, wt)

wt = [1,50,1]
map_square(plt, 3, A, wt)
map_square(plt, 2, A, wt)
map_square(plt, 1, A, wt)

wt = [1,1,50]
map_square(plt, 3, A, wt)
map_square(plt, 2, A, wt)
map_square(plt, 1, A, wt)

##############
## Problem 1
##############

## need to implement the projective moment map

## d = 2, n = 3
A = hcat([[0, 0], [0, 1], [1, 0]]...)
A′ = vcat(A, [1 1 1])
T = convex_hull(transpose(A)) ## triangle
X = toric_ideal(transpose(A′)) ## == ⟨ 0 ⟩ because X = ℙ²
dim(X) == 3

wt = [1,1,1]
ϕ = param(A)
ψ = proj_mm(A, wt)
f = ψ ∘ ϕ

plt = plot(xlims=[0,1], ylims=[0,1])
α = 100
pts = sample_square(α, 2)
image_pts = hcat(f.(pts)...)
scatter!(plt, image_pts[1,:], image_pts[2,:], label="$(α)")

##############
## Problem 2
##############

A = hcat([[1, 0], [1, 1], [1, 2]]...)
A = hcat([[1, -1], [1, 0], [1,1]]...) ## does this make the parametrization better
IY = toric_ideal(transpose(A))
ϕ = param(A)

plt = plot()
pts = hcat([ϕ(p) for p in sample_square(-10, 10, 2)]...)
scatter!(plt, pts[1,:], pts[2,:], pts[3,:], label="toric variety")
pos_pts = hcat([ϕ(p) for p in sample_square(10, 2)]...)
scatter!(plt, pos_pts[1,:], pos_pts[2,:], pos_pts[3,:], camera = (10, 20), label="positive points")