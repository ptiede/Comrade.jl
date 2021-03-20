using Memoization
using StatsBase
using LRUCache
const rad2μas = 180.0/π*3600*1e6

function nuts_sample(start, lj, forward_lj, nsample, nwarmup)
    metric = DiagEuclideanMetric(length(start))
    hamiltonian = Hamiltonian(metric, lj, forward_lj)
    prop = AdvancedHMC.NUTS(Leapfrog(find_good_stepsize(hamiltonian, start)))
    initial_step = find_good_stepsize(hamiltonian, start)
    integrator = Leapfrog(initial_step)
    adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator))
    samples, stats = sample(hamiltonian,
                            prop,
                            start, nsample,
                            adaptor, nwarmup;
                            progress=true, drop_warmup=true)
    return samples, stats
end


function threaded_nlopt(nopt, fopt, srange,  maxevals; starts=nothing)
    results = [zeros(lj.transform.dimension) for _ in 1:nopt]
    divs = zeros(nopt)
    lower = first.(srange)
    upper = last.(srange)
    #starts = rand(truncated(Normal(0.0, 3.0), first(srange[1]), last(srange[1])),
    #              length(srange), nopt)
    Threads.@threads for i in 1:nopt
        nlopt = NLopt.Opt(:LD_LBFGS, length(srange))
        lower_bounds!(nlopt, lower)
        upper_bounds!(nlopt, upper)
        max_objective!(nlopt, fopt)
        xtol_rel!(nlopt,1e-12)
        maxeval!(nlopt, maxevals)
        start = rand(length(srange)).*(upper - lower) + lower
        (minf,minx,ret) = NLopt.optimize(nlopt, start)
        results[i] = minx
        divs[i] = minf
    end
    I = sortperm(divs; rev=true)

    return results[I], divs[I]
end


function threaded_bbopt(nopt, lj, srange, maxevals)
    results = [zeros(lj.transform.dimension) for _ in 1:nopt]
    divs = zeros(nopt)

    Threads.@threads for i in 1:nopt
        res = bboptimize(x->-lj(x), SearchRange=srange, MaxFuncEvals=maxevals, TraceMode=:compact)
        results[i] = best_candidate(res)
        divs[i] = best_fitness(res)
    end
    I = sortperm(divs)
    return results[I], divs[I]
end

function plot_mean(sims, imi)
    p = plot()
    mean_img = mean(getproperty.(sims, :im))
    println(size(mean_img))
    heatmap!(p,
              pixel_iterator(sims[1])...,
              mean_img,
              aspect_ratio=:equal,
              size=(500,400),
              xflip=true
            )
    xlims!(p,-80.0, 80.0)
    ylims!(p,-80.0, 80.0)
    title!(p,"Mean Image Frame: $imi")
    return p
end

function plot_samples(sims, rchi2, imi)
    anim = @animate for s in sims[1:end]
        heatmap(pixel_iterator(s)..., s.im, aspect_ratio=:equal, size=(500,400), xflip=true)
        xlims!(-80.0,80.0)
        ylims!(-80.0, 80.0)
        title!("Frame $imi,     χ²ᵣ = $rchi2")
    end
    return anim
end


function threaded_opt(nopt, lj, srange, maxevals)
    results = [zeros(lj.transform.dimension) for _ in 1:nopt]
    divs = zeros(nopt)
    lower = first.(srange)
    upper = last.(srange)
    Threads.@threads for i in 1:nopt
        try
            start = lower .+ (upper .- lower).*rand(lj.transform.dimension)
            res = optimize(x->-lj(x), start, LBFGS(), autodiff=:forward)
            results[i] = Optim.minimizer(res)
            divs[i] = Optim.minimum(res)
        catch
            divs[i] = Inf
        end
    end
    I = sortperm(divs)
    return results[I], divs[I]
end

function chi2(m::ROSE.AbstractModel, data; ferr2 = 0.0, floor = 0.0)
    u = data.u[:,1]
    v = data.v[:,1]
    vdata = StructArray{Complex{Float64}}(re=data.visr[:,1], im=data.visi[:,1])
    visi = data.visi
    vmodel = ROSE.visibilities(m, u/rad2μas, v/rad2μas)
    err = sqrt.(data.error[:,1].^2 + ferr2*abs2.(vmodel) .+ floor^2)
    nres = (vmodel-vdata)./err
    chi2 = sum(abs2, (vmodel-vdata)./err)
    return chi2, nres
end

function split_data(obs, dt)
    starttime = getdata(obs, :time)[1]
    endtime = getdata(obs, :time)[end]
    sdata = []
    times = collect(range(starttime, stop=endtime, step=dt/3600))
    obstimes = getdata(obs,:time)
    for i in 2:length(times)
        tmp = obs.data[times[i-1].<= obstimes .< times[i]]
        if length(tmp) > 30
            push!(sdata, tmp)
        end
    end
    return sdata
end



function make_sims(chain, range, nx, ny)
    subchain = chain[range]
    coeffs = getproperty.(subchain, :coeffs)
    ϵ = getproperty.(subchain, :ϵ)
    scx = getproperty.(subchain, :scx)
    scy = getproperty.(subchain, :scy)
    f = getproperty.(subchain, :f)
    x0 = getproperty.(subchain, :x0)
    y0 = getproperty.(subchain, :y0)
    ϵ = getproperty.(subchain, :ϵ)
    rimg = @. shifted.(renormed.(stretched.(RImage.(reshape.(coeffs, nx,ny), SqExpKernel.(ϵ)),
                     scx, scy
                    ), f),
                    x0, y0
                )
    cres = getproperty.(subchain, :ring)
    router = getproperty.(cres, :R)
    rinner = router.*(1.0 .- getproperty.(cres, :ψ))
    tmp = ConcordanceCrescent.(router, rinner, 0.0, getproperty.(cres, :s))
    cres = shifted.(
                ROSE.renormed.(
                    #smoothed.(
                    rotated.(
                        tmp, getproperty.(cres,:ξ)),
                        #getproperty.(cres, :σ)),
                    getproperty.(cres,:f)),
                getproperty.(cres,:x0), getproperty.(cres,:y0)
                )

    ims = rimg .+ cres
    sims = ROSE.stokesimage.(ims, 512, 512, 160.0, 160.0)
end


function make_static_sims(chain, range, nx, ny; scx=70.0, scy=70.0)
    subchain = chain[range]
    coeffs = getproperty.(subchain, :coeffs)
    f = getproperty.(subchain, :f)
    ϵ = getproperty.(subchain, :ϵ)
    ims = renormed.(stretched.(RImage.(reshape.(coeffs, nx,ny), SqExpKernel.(ϵ)),
                     scx, scy), f)

    sims = ROSE.stokesimage.(ims, 128, 128, 160.0, 160.0)
    return sims
end


@memoize function reversecache(f, dim)
    ftape = ReverseDiff.GradientTape(f, rand(dim))
    return ReverseDiff.compile(ftape)
end

@memoize LRU(maxsize=1) function forwardchunk(f, dim, csize)
    cfg = ForwardDiff.GradientConfig(f, rand(dim), ForwardDiff.Chunk{csize}())
    return cfg
end

function make_geometric_sims(chain, range, model)
    subchain = chain[range]

    router = getproperty.(subchain, :R)
    rinner = router .* (1 .-getproperty.(subchain,:ψ))
    shift = getproperty.(subchain,:τ).*(router .- rinner)
    img = ConcordanceCrescent.(router, rinner, shift, getproperty.(subchain,:s))
    mopt = shifted.(
            renormed.(
                rotated.(
                    ROSE.smoothed.(img, getproperty.(subchain,:σ)),
                getproperty.(subchain,:ξ)),
            getproperty.(subchain,:f)),
        getproperty.(subchain,:x0), getproperty.(subchain,:y0))


    # Gaussian parameters
    if model === cres1g || model === cres2g
        g1chain = getproperty.(subchain, :g1)
        scx = getproperty.(g1chain,:σg)./sqrt.(1 .- getproperty.(g1chain,:τg))
        scy = getproperty.(g1chain,:σg)./sqrt.(1 .+ getproperty.(g1chain,:τg))
        mg1 = shifted.(renormed.(rotated.(stretched.(Ref(Gaussian()), scx, scy),
                    getproperty.(g1chain,:ξg) ), getproperty.(g1chain,:fg)),
                    getproperty.(g1chain,:xg), getproperty.(g1chain,:yg))
        mopt = mopt .+ mg1
    end

    if model === cres2g
        g2chain = getproperty.(subchain, :g2)
        scx = getproperty.(g2chain,:σg)./sqrt.(1 .- getproperty.(g2chain,:τg))
        scy = getproperty.(g2chain,:σg)./sqrt.(1 .+ getproperty.(g2chain,:τg))
        mg2 = shifted.(renormed.(rotated.(stretched.(Ref(Gaussian()), scx, scy),
                        getproperty.(g2chain,:ξg) ), getproperty.(g2chain,:fg)),
                        getproperty.(g2chain,:xg), getproperty.(g2chain,:yg))
        mopt = mopt .+ mg2
    end
    sims = ROSE.stokesimage.(mopt, 256, 256, 240.0, 240.0)
    println("length of sims: $(length(sims))")
    return sims
end


function plot_vis_comp(data, vmodel; ferr2=0.0, floor=0.0)
    p = plot()
    error = sqrt.(data.error.^2 + ferr2*abs2.(vmodel) .+ floor^2)
    scatter!(p,hypot.(data.u, data.v)/1e9, data.visr,
                  yerr=error,
                  color=:cornflowerblue, label="Data Real", markershape=:square, alpha=0.5)
    scatter!(p,hypot.(data.u, data.v)/1e9, real.(vmodel),
                  color=:blue, label="Model Real")
    scatter!(p,hypot.(data.u, data.v)/1e9, data.visi,
                  yerr=error,
                  color=:orange, label="Data Imag", markershape=:square, alpha=0.5)
    scatter!(p,hypot.(data.u, data.v)/1e9, imag.(vmodel),
                  color=:red, label="Model Imag")
    xlabel!(p,"uv dist (Gλ)")
    ylabel!(p,"Complex Vis (Jy)")
    return p
end


function plot_res_uv(data, nres, rchi2, imi)
    p = plot()
    scatter!(p,hypot.(data.u,data.v)/1e9, real.(nres), label="Real Vis")
    scatter!(p,hypot.(data.u,data.v)./1e9, imag.(nres), label="Imag Vis")
    xlabel!(p,"uv dist (Gλ)")
    ylabel!(p,"Normalized Residuals")
    ylims!(p, -6.0, 6.0)
    title!(p, "Frame $imi,     χ²ᵣ = $rchi2")
    return p
end

function plot_res_time(data, nres, rchi2, imi)
    p = plot()
    scatter!(p,data.time, real.(nres), label="Real Vis")
    scatter!(p,data.time, imag.(nres), label="Imag Vis")
    xlabel!(p,"time [hr]")
    ylabel!(p,"χᵣ²")
    ylims!(p,-6.0, 6.0)
    title!(p,"Frame $imi,     χ²ᵣ = $rchi2")
    return p
end

function plot_res_density(nres)
    p = plot()
    density!(p,real.(nres), label="Real")
    density!(p,imag.(nres), label="Imag")
    plot!(p,x->pdf(Normal(), x), label="Std. Normal")
    xlabel!(p,"Normalized Residuals")
    return p
end


struct LogJoint{D,M,T}
    data::D
    model::M
    transform::T
end

trans(lj::LogJoint{D,M,T}) where {D,M,T} = lj.transform
dim(lj) = trans(lj).dimension


function LogJoint(data, model)
    return LogJoint(data, model, xform(model,data))
end

function (ℓ::LogJoint)(x)
    (θ, logjac) = Soss.transform_and_logjac(trans(ℓ), x)
    return logpdf(ℓ.model, merge(θ, ℓ.data)) + logjac
end


"""
    Importance samples the nested sampling chain
    with weights defined in the chain object.

    While nsamp can be whatever you want, please be careful!
    There are only a finte number of effective samples in the chain.
    Going much farther than the ESS doesn't lower the MCSE.
"""
function importance_sample(chain, nsamp::Int)
    weights = ProbabilityWeights(chain[:weights].data[:,1])
    return sample(chain[:,1:end-1,:], weights, nsamp)
end

"""
    ns_ess(chain)
Computes the effective sample size (ESS) of
the nested sampling run. This isn't the best
estimate of the ess but it will give you reasonable
results.
"""
function ns_ess(chain)
    return sum(chain[:weights])^2/sum(chain[:weights].^2)
end
