using Pkg; Pkg.activate(@__DIR__)
using Comrade
using Distributions
using ComradeOptimization
using ComradeAHMC
using GalacticBBO
using Plots
using StatsBase
using GalacticOptimJL
using DistributionsAD

# load eht-imaging we use this to load eht data
load_ehtim()
# To download the data visit https://doi.org/10.25739/g85n-f134
obs = ehtim.obsdata.load_uvfits(joinpath(@__DIR__, "SR1_M87_2017_096_hi_hops_netcal_StokesI.uvfits"))
obs.add_scans()
# kill 0-baselines since we don't care about
# large scale flux and make scan-average data
obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.02)
# extract log closure amplitudes and closure phases
dlcamp = extract_lcamp(obs)
dcphase = extract_cphase(obs)
lklhd = RadioLikelihood(dlcamp, dcphase)

# Build the Model. Here we we a struct to hold some caches
# which will speed up imaging
struct ImModel{C}
    cache::C
    fov::Float64
    npix::Int
end

# For our model we will be using a rasterized image. This can be viewed as something like a
# non-parametric model. As a result of this we will need to use a `modelimage` object to
# store cache information we will need to compute the numerical FT.
function (model::ImModel)(θ)
    (;c) = θ
    #Construct the image model
    cmat = reshape(c, model.npix, model.npix)
    img = IntensityMap(cmat, model.fov, model.fov, BSplinePulse{3}())
    #Create the modelimage object that will use a cache to compute the DFT
    m = modelimage(img, cache)
end

npix = 32
fovxy = μas2rad(65.0)
# Now we can feed in the array information to form the cache. We will be using a DFT since
# it is efficient for so few pixels
cache = create_cache(Comrade.DFTAlg(dlcamp), IntensityMap(rand(npix,npix), fovxy, fovxy, BSplinePulse{3}()))
mms = ImModel(cache, fovxy, npix)
# We will use a Dirichlet prior to enforce that the flux sums to unity since closures are
# degenerate to total flux.
prior = (c = MvNormal(fill(-5.0, npix^2), 1.0),)

post = Posterior(lklhd, prior, mms)
tpost = asflat(post)

# Let's run an optimizer to get a nice starting location
# It turns out that gradients are really helpful here
ndim = dimension(tpost)
f = OptimizationFunction(tpost, GalacticOptim.AutoZygote())
x0 = rand(ndim)
prob = OptimizationProblem(f, x0, nothing)
sol = solve(prob, LBFGS(); maxiters=4_000, show_trace=true, show_every=50, g_tol=1e-2)

xopt = transform(tpost, sol)

# Let's see how the fit looks
residual(mms(xopt), dlcamp)
residual(mms(xopt), dcphase)
plot(mms(xopt), fovx=1.2*fovxy, fovy=1.2*fovxy)


# now we sample using hmc
metric = DenseEuclideanMetric(ndim)
hchain, stats = sample(post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 3000; nadapts=2000, init_params=xopt)

lca = logdensityof(lklhd.lklhds[1])
lcp = logdensityof(lklhd.lklhds[2])

foox = let lca=lca, lcp=lcp, pr=tpost.lpost.prior, tr=tpost.transform, plan = cache.plan, phases=cache.phases, dmat=dlcamp.config.designmat, dmatc=dcphase.config.designmat
    x->begin
        y = transform(tr, x)
        lp = logdensityof(pr, y)
        vis = plan*exp.(y.c).*phases
        return lca(vis) + lcp(vis) + lp
    end
end




function bench(npix, lklhd, dlcamp, )
    fovxy = μas2rad(65.0)
    # Now we can feed in the array information to form the cache. We will be using a DFT since
    # it is efficient for so few pixels
    cache = create_cache(Comrade.DFTAlg(dlcamp), IntensityMap(rand(npix,npix), fovxy, fovxy, BSplinePulse{3}()))
    mms = ImModel(cache, fovxy, npix)
    # We will use a Dirichlet prior to enforce that the flux sums to unity since closures are
    # degenerate to total flux.
    prior = (c = MvNormal(fill(-5.0, npix^2), 1.0),)

    post = Posterior(lklhd, prior, mms)
    tpost = asflat(post)

    # Let's run an optimizer to get a nice starting location
    # It turns out that gradients are really helpful here
    ndim = dimension(tpost)

    lca = logdensityof(lklhd.lklhds[1])
    lcp = logdensityof(lklhd.lklhds[2])


    foox = let lca=lca, lcp=lcp, pr=tpost.lpost.prior, tr=tpost.transform, plan = cache.plan, phases=cache.phases, dmat=dlcamp.config.designmat, dmatc=dcphase.config.designmat
        x->begin
            y = transform(tr, x)
            lp = logdensityof(pr, y)
            vis = plan*exp.(y.c).*phases
            return lca(vis) + lcp(vis) + lp
        end
    end
    x0 = randn(ndim)
    ftime = @belapsed $(foox)($x0)
    gftime = @belapsed Zygote.gradient($foox, $x0)
    gfftime = @belapsed AD.gradient($(AD.ForwardDiffBackend{npix}()), $foox, $x0)
    return npix, ftime, gftime, gftime/ftime, gfftime, gfftime/ftime
end


btimes = [bench(npix, lklhd, dlcamp) for npix in 4:2:32]

using Plots

f = scatter(getindex.(btimes, 1).^2, getindex.(btimes, 4), yscale=:log10, label="Reverse Mode")
scatter!(f, getindex.(btimes, 1).^2, getindex.(btimes, 6), label="Forward Mode")

ylabel!(f, "∇f/f")
xlabel!(f, "npix²")
xlims!(4, 512)
ylims!(1.0, 10000.0)
savefig(f, "image_scaling_grad_order.png")


ndims = getindex.(btimes, 1).^2
f2 = scatter(ndims, getindex.(btimes, 2), yscale=:log10, label="F(x)")
scatter!(f2, ndims, getindex.(btimes, 3), label="Reverse ∇F(x)")
scatter!(f2, ndims, getindex.(btimes, 5), label="Forward ∇F(x)")
plot!(f2, ndims, ndims*btimes[1][2]/200 .+ ndims[1]*btimes[1][2]/20, color=:blue, label="O(npix²)")
plot!(f2, ndims, ndims*btimes[1][2]/200*10 .+ ndims[1]*btimes[1][2]/20*10, color=:orange, label="10⋅O(npix²)")
plot!(f2, ndims, ndims.*ndims*btimes[1][2]/50 .+ ndims[1]*ndims[1]*btimes[1][2]/10, color=:green, label="O(npix⁴)")
ylabel!(f2, "Runtime (s)")
xlabel!(f2, "npix²")
xlims!(f2, 4, 512)
savefig(f2, "image_scaling_grad_times.png")

function logistic_logjac(x)
    mx = -abs(x)
    mx - 2*log1p(exp(mx))
end



function simplex_fwd(y)
    x = similar(y, length(y)+1)
    logjac = zero(eltype(y))
    stick = one(eltype(y))
    n = length(y)+1
    @inbounds for i in eachindex(y)
        z = logistic(y[i] + log(n-i))
        x[i] = stick*(1-z)
        stick *= z
        logjac += log(abs(x[i]/(z-1))) - logistic_logjac(z)
    end
    x[end] = stick

    return x, logjac
end

using ChainRulesCore

function ChainRulesCore.rrule(::typeof(simplex_fwd), y)
    x, ℓ = simplex_fwd(y)
    function _simplex_fwd_pullback(ΔX)
        f̄ = NoTangent()
        (Δx, Δℓ) = ΔX
        jx = similar(x, length(x), length(y))
        for i in eachindex(y)

        end
    end
end

function mymul!(R, A, B)
    @assert axes(A,2) == axes(B,1)
    @inbounds @simd for i in eachindex(R)
        R[i] = 0
    end
    @inbounds for j in axes(B, 2), i in axes(A, 1)
        @inbounds @simd for k in axes(A,2)
            R[i,j] += A[i,k] * B[k,j]
        end
    end
    nothing
end

A = rand(5, 3)
B = rand(3, 7)

R = zeros(size(A,1), size(B,2))
∂z_∂R = rand(size(R)...)  # Some gradient/tangent passed to us

∂z_∂A = zero(A)
∂z_∂B = zero(B)

Enzyme.autodiff(mymul!, Const, Duplicated(R, ∂z_∂R), Duplicated(A, ∂z_∂A), Duplicated(B, ∂z_∂B))



y0 = randn(16^2)
∂y0 = zero(y0)

x0 = zeros(length(y0)+2)
∂x0 = zeros(size(x0)...)

d = Dirichlet(npix^2, 1.0)


fooE = let lca=dlcamp, cp=dcphase, pr=d, plan = cache.plan, phases=cache.phases, dmat=dlcamp.config.designmat, dmatc=dcphase.config.designmat
    x->begin
        y, lj = simplex_fwd(x)
        lp = logdensityof(pr, y)
        vis = plan*y.*phases
        mlca = dmat*log.(abs.(vis))
        mcp = dmatc*angle.(vis)
        l1 = sum(abs2, (lca[:amp] .- mlca)./lca[:error])
        l2 = sum(abs2, (cp[:phase] .- mcp)./cp[:error])
        return -0.5*(l1 + l2) + lp + lj
    end
end


Enzyme.autodiff(simplex_fwd!, Const, Duplicated(x0, ∂x0), Duplicated(y0, ∂y0))

ℓ = logdensityof(d)

struct LT{D}
    ℓ::D
end

function (lt::LT)(z)
    x, logjac = simplex_fwd(z)
    return lt.ℓ(x) + logjac
end

ℓt = LT(ℓ)

∂z0 = zero(z0)
@benchmark Enzyme.gradient($(Enzyme.ReverseMode()), $ℓt, $z0)
@benchmark Zygote.gradient(ℓt, z0)

function mymul!(R, A, B)
    @assert axes(A,2) == axes(B,1)
    @inbounds @simd for i in eachindex(R)
        R[i] = 0
    end
    @inbounds for j in axes(B, 2), i in axes(A, 1)
        @inbounds @simd for k in axes(A,2)
            R[i,j] += A[i,k] * B[k,j]
        end
    end
    nothing
end


A = rand(5, 3)
B = rand(3, 7)

R = zeros(size(A,1), size(B,2))
∂z_∂R = rand(size(R)...)  # Some gradient/tangent passed to us

∂z_∂A = zero(A)
∂z_∂B = zero(B)

Enzyme.autodiff(mymul!, Const, Duplicated(R, ∂z_∂R), Duplicated(A, ∂z_∂A), Duplicated(B, ∂z_∂B))
