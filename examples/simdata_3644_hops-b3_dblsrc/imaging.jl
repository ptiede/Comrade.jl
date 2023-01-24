#LinearAlgebra.BLAS.set_number_threads(20)
using Pkg;
Pkg.add(url="https://github.com/ptiede/RadioImagePriors.jl")
using Comrade
using Distributions
using DistributionsAD
using ComradeOptimization
using ComradeAHMC
using OptimizationBBO
using Zygote
using Plots
using StatsBase
using OptimizationOptimJL
using RadioImagePriors
using TupleVectors

using Glob 
uvf_list= glob("*.uvfits")
dir_list = String[]
for s in uvf_list
	x = replace(s, r".uvfits$"=>"")
	mkdir(x)
	push!(dir_list,x)

	# load eht-imaging we use this to load eht data
	load_ehtim()
	obs = ehtim.obsdata.load_uvfits(s)
	obs.add_scans()
	# make scan-average data #.flag_uvdist(uv_min=0.1e9)
	obs = obs.flag_uvdist(uv_min=0.1e9).avg_coherent(0.0, scan_avg=true).add_fractional_noise(0.01)

	# extract amplitudes and closure phases
	damp = extract_amp(obs)
	dcphase = extract_cphase(obs)
	lklhd = RadioLikelihood(damp, dcphase)

	# Build the Model. Here we make a struct to hold some caches
	# This will be useful to hold precomputed caches

	struct GModel{C,G}
		cache::C
		gcache::G
		fovx::Float64
		fovy::Float64
		npixx::Int
		npixy::Int
		function GModel(obs::Comrade.EHTObservation, fovx, fovy, npixx, npixy)
		    buffer = IntensityMap(zeros(npixx, npixy), fovx, fovy, BSplinePulse{3}())
		    cache = create_cache(DFTAlg(obs), buffer)
		    gcache = GainCache(scantable(obs))
		    return new{typeof(cache), typeof(gcache)}(cache, gcache, fovx, fovy, npixx, npixy)
		end
	end
	
	function (model::GModel)(θ)
		(;c, f, lgamp) = θ
		# Construct the image model
		img = IntensityMap(f*c, model.fovx, model.fovy, BSplinePulse{3}())
		m = modelimage(img, model.cache)
		# Now corrupt the model with Gains
		g = exp.(lgamp)
		Comrade.GainModel(model.gcache, g, m)
	end
		

	# First we define the station gain priors
	distamp = (AA = Normal(0.0, 0.1),
		       AX = Normal(0.0, 0.1),
		       LM = Normal(0.0, 0.9),
		       SW = Normal(0.0, 0.1),
		       MM = Normal(0.0, 0.1),
		       PV = Normal(0.0, 0.1),
		       MG = Normal(0.0, 0.1),
		       GL = Normal(0.0, 0.5)       
		       )
	
	npixx = 13
	npixy = 11
	fovx = μas2rad(65.0)
	fovy = μas2rad(55.0)
	
	prior = (
		      c = ImageDirichlet(0.5, npixx, npixy),
		      f = Uniform(0.4, 0.6),
		      lgamp = Comrade.GainPrior(distamp, scantable(damp)),
		    )


	mms = GModel(damp, fovx, fovy, npixx, npixy)	

	post = Posterior(lklhd, prior, mms)

	tpost = asflat(post)

	# We will use HMC to sample the posterior.
	# First to get in the right ballpark we will use `BlackBoxOptim.jl`
	
	############################
	
	#ndim = dimension(tpost)
	#f = OptimizationFunction(tpost, Optimization.AutoZygote())
	#prob = OptimizationProblem(f, rand(ndim), nothing, lb=fill(-5.0, ndim), ub=fill(5.0, ndim))
	#sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000_000)
	############################
	
	# Now lets zoom to the peak using LBFGS
	ndim = dimension(tpost)
	using Zygote
	f = OptimizationFunction(tpost, Optimization.AutoZygote())
	prob = OptimizationProblem(f, randn(ndim), nothing)
	ℓ = logdensityof(tpost)
	sol = solve(prob, LBFGS(); maxiters=8_000, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)

	xopt = transform(tpost, sol)
	
	p4 = plot(mms(xopt), fovx=fovx, fovy=fovy, dpi=500)
	display(p4)
	savefig("./"*x*"/map.png")
	
	# now we sample using hmc
	metric = DiagEuclideanMetric(ndim)
	hchain, stats = sample(post, AHMC(;metric, autodiff=AD.ZygoteBackend()), 4000; nadapts=3000, init_params=xopt)

	using Serialization
	serialize("./"*x*"/"*x*"_fov_65x55_npix_13x11.jls",
		           Dict(:chain=>hchain,
		               :stats=>stats,
		               :xopt=>xopt,
		               :npixx => npixx,
		               :npixy => npixy,
		               :fovx => fovx,
		               :fovy => fovy
		           ))

	# Plot the mean image and standard deviation image
	using StatsBase
	samples = mms.(sample(hchain, 500))
	imgs = intensitymap.(samples, fovx, fovy, 96, 96)

	mimg, simg = mean_and_std(imgs)

	p1 = plot(mimg, title="Mean", clims=(0.0, maximum(mimg)), dpi=500)
	display(p1)
	savefig("./"*x*"/mean.png")
	p2 = plot(simg,  title="Std. Dev.", clims=(0.0, maximum(mimg)), dpi=500)
	display(p2)
	savefig("./"*x*"/std.png")
	p3 = plot(simg./mimg,  title="Fractional Error", xlims=(-32.5,32.5), ylims=(-27.5,27.5), dpi=500)
	display(p3)
	savefig("./"*x*"/error.png")

	p5 = residual(mms(xopt), dcphase, ylabel="Closure Phase Norm. Res.", dpi=500)
	display(p5)
	savefig("./"*x*"/res_cphase.png")
	p6 = residual(mms(xopt), damp, ylabel="Amp. Norm. Res.", dpi=500)
	display(p6)
	savefig("./"*x*"/res_amp.png")

	#Plot the gain table with error bars
	gamps = exp.(hcat(hchain.lgamp...))
	mga = mean(gamps, dims=2)
	sga = std(gamps, dims=2)

	using Measurements
	gmeas = measurement.(mga, sga)
	ctable = caltable(mms.gcache, vec(gmeas))
	p7 = plot(ctable, layout=(3,3), size=(600, 500), datagains=true, dpi=500)
	display(p7)
	savefig("./"*x*"/gains.png")
	
	p8 = plot(mms(xopt), dcphase, dpi=500)
	display(p8)
	savefig("./"*x*"/cphase.png")
	
	p9 = plot(mms(xopt), dcphase, xlims=(0,1.5e8), dpi=500)
	display(p9)
	savefig("./"*x*"/cphase2.png")
end

