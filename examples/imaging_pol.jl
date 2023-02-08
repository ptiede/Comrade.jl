# # Polarized Image and Instrumental Modeling

# In this tutorial we will analyze an simulated simple polarized dataset to demonstrate
# Comrade's polarized imaging capabilities.

# ## Introduction to Polarized Imaging


using Pkg; Pkg.activate(@__DIR__)

using Comrade

# ## Load the Data
# To download the data visit https://doi.org/10.25739/g85n-f134
# To load the eht-imaging obsdata object we do:
obs = load_ehtim_uvfits(joinpath(@__DIR__, "PolarizedExamples/polarized_gaussian.uvfits"),
                        joinpath(@__DIR__, "PolarizedExamples/array.txt"))

# Now we do some minor preprocessing:
#   - Scan average the data since the data have been preprocessed so that the gain phases
#      coherent.
#   - Add 1% systematic noise to deal with calibration issues that cause 1% non-closing errors.
obs = scan_average(obs)

# Now we extract our coherency matrices
dvis = extract_coherency(obs)

# ##Building the Model/Posterior

# Now we must build our intensity/visibility model. That is, the model that takes in a
# named tuple of parameters and perhaps some metadata required to construct the model.
# For our model we will be using a raster or `ContinuousImage` for our image model.
# Unlike other imaging examples
# (e.g., [Imaging a Black Hole using only Closure Quantities](@ref)) we also need to include
# a model for the intrument, i.e., gains as well. The gains will be broken into two components
#   - Gain amplitudes which are typically known to 10-20% except for LMT which has large issues
#   - Gain phases which are more difficult to constrain and can shift rapidly.
# The model is given below:


function model(θ, metadata)
    (;c, f, p, angparams, dRx, dRy, dLx, dLy, lgp, gpp, lgr, gpr) = θ
    (; grid, cache, pulse, tcache, gcache, dcache) = metadata
    # Construct the image model
    # produce Stokes images from parameters
    imgI = f*c
    # Converts from poincare sphere parameterization of polzarization to Stokes Parameters
    pimg = PoincareSphere2Map(imgI, p, angparams, grid)
    cimg = ContinuousImage(pimg, pulse)

    m = modelimage(cimg, cache)
    jT = jonesT(tcache)
    # calibration parameters
    # Gain product parameters
    gP = cis.(lgp./2 .+ 1im.*gpp./2)
    Gp = jonesG(gP, gP, gcache)
    # Gain ratio
    gRat = cis.(lgr./2 .+ 1im.*gpr./2)
    Gr = jonesG(gRat, inv.(gRat), dcache)
    D = jonesD(complex.(dRx, dRy), complex.(dLx, dLy), dcache)
    J = Gp*Gr*D*jT
    return JonesModel(J, m, tcache)
end





# First we define the station gain priors
distamp = (AA = Normal(0.0, 0.1),
           AP = Normal(0.0, 0.1),
           LM = Normal(0.0, 0.3),
           AZ = Normal(0.0, 0.1),
           JC = Normal(0.0, 0.1),
           PV = Normal(0.0, 0.1),
           SM = Normal(0.0, 0.1),
           )


distphase = (AA = DiagonalVonMises([0.0], [inv(1e-4)]),
             AP = DiagonalVonMises([0.0], [inv(π^2)]),
             LM = DiagonalVonMises([0.0], [inv(π^2)]),
             AZ = DiagonalVonMises([0.0], [inv(π^2)]),
             JC = DiagonalVonMises([0.0], [inv(π^2)]),
             PV = DiagonalVonMises([0.0], [inv(π^2)]),
             SM = DiagonalVonMises([0.0], [inv(π^2)]),
           )

distD = ( AA = Normal(0.0, 0.1),
          AP = Normal(0.0, 0.1),
          LM = Normal(0.0, 0.1),
          AZ = Normal(0.0, 0.1),
          JC = Normal(0.0, 0.1),
          PV = Normal(0.0, 0.1),
          SM = Normal(0.0, 0.1),
        )


# Set up the cache structure
fovx = μas2rad(30.0)
fovy = μas2rad(30.0)
nx = 3
ny = floor(Int, fovy/fovx*nx)

grid = imagepixels(fovx, fovy, nx, ny)
buffer = IntensityMap(zeros(nx, ny), grid)
pulse = BSplinePulse{3}()
cache = create_cache(NFFTAlg(dvis), buffer, pulse)
tcache = TransformCache(dvis; add_fr=false)
gcache = JonesCache(dvis, ScanSeg())
dcache = JonesCache(dvis, TrackSeg())
metadata = (;cache, grid, pulse, tcache, gcache, dcache)


prior = (
          c = ImageDirichlet(1.0, nx, ny),
          f = Uniform(0.7, 1.3),
          p = ImageUniform(nx, ny),
          angparams = ImageSphericalUniform(nx, ny),
          lgp = CalPrior(distamp, gcache),
          gpp = CalPrior(distphase, gcache),
          lgr = CalPrior(distamp, dcache),
          gpr = CalPrior(distphase, dcache),
          dRx = CalPrior(distD, dcache),
          dRy = CalPrior(distD, dcache),
          dLx = CalPrior(distD, dcache),
          dLy = CalPrior(distD, dcache),
          )



lklhd = RadioLikelihood(model, metadata, dvis)

post = Posterior(lklhd, prior)

tpost = asflat(post)
ndim = dimension(tpost)

ℓ = logdensityof(tpost)



using Zygote
f = OptimizationFunction(tpost, Optimization.AutoZygote())
prob = OptimizationProblem(f, prior_sample(tpost), nothing)
sol = solve(prob, LBFGS(); maxiters=8_000, callback=(x,p)->(@info ℓ(x); false), g_tol=1e-1)
xopt = transform(tpost, sol)


# Let's see how the fit look. First let's load the ground truth image to compare everything
imgtrue = Comrade.load(joinpath(@__DIR__, "PolarizedExamples/polarized_gaussian.fits"), StokesIntensityMap)

# Now let's see what our MAP image looks like
img = intensitymap(model(xopt, metadata), fovxy, fovxy, npix*5, npix*5)
plot(img)
Comrade.save(joinpath(@__DIR__, "test.fits"), img)

residual(model(xopt, metadata), dvis)
plot(model(xopt, metadata), dvis)

# Let's also plot the calibration tables for gains
gL = Comrade.caltable(gcache, exp.(xopt.lgp .+ xopt.lgr))
plot(gL, layout=(3,3), size=(600,500))

gR = Comrade.caltable(gcache, exp.(xopt.lgR))
plot!(gR, layout=(3,3), size=(600,500))

# Let's also plot the calibration tables for gains
gpL = Comrade.caltable(gcache, (xopt.gpL))
plot(gpL, layout=(3,3), size=(600,500))

gpR = Comrade.caltable(gcache, (xopt.gpR))
plot!(gpR, layout=(3,3), size=(600,500))

# And the calibration tables for d-terms
dR = caltable(dcache, complex.(xopt.dRx, xopt.dRy))
dL = caltable(dcache, complex.(xopt.dLx, xopt.dLy))



vis = dvis[:measurement]
err = dvis[:error]
mvis = visibilities(model(xopt, metadata), arrayconfig(dvis))
uvdist = hypot.(values(getuv(dvis))...)

p1 = scatter(uvdist, real.(getindex.(vis, 1, 1)), yerr=getindex.(err, 1, 1), title="RR", color=:blue, marker=:circ, label="Data real")
scatter!(uvdist, imag.(getindex.(vis, 1, 1)), yerr=getindex.(err, 1, 1), color=:orange, marker=:circ, label="Data imag")
scatter!(uvdist, real.(getindex.(vis, 1, 1)), color=:cyan, alpha=0.5, marker=:square, label="Model real")
scatter!(uvdist, imag.(getindex.(vis, 1, 1)), color=:yellow, alpha=0.5, marker=:square, label="Model imag")

p2 = scatter(uvdist, real.(getindex.(vis, 2, 1)), yerr=getindex.(err, 2, 1), title="LR", color=:blue, marker=:circ, label="Data real RR")
scatter!(uvdist, imag.(getindex.(vis, 2, 1)), yerr=getindex.(err, 2, 1), color=:orange, marker=:circ, label="Data imag RR")
scatter!(uvdist, real.(getindex.(vis, 2, 1)), color=:cyan, alpha=0.5, marker=:square, label="Model real RR")
scatter!(uvdist, imag.(getindex.(vis, 2, 1)), color=:yellow, alpha=0.5, marker=:square, label="Model imag RR", legend=:false)

p3 = scatter(uvdist, real.(getindex.(vis, 1, 2)), yerr=getindex.(err, 1, 2), title="RL", color=:blue, marker=:circ, label="Data real RR")
scatter!(uvdist, imag.(getindex.(vis, 1, 2)), yerr=getindex.(err, 1, 2), color=:orange, marker=:circ, label="Data imag RR")
scatter!(uvdist, real.(getindex.(vis, 1, 2)), color=:cyan, alpha=0.5, marker=:square, label="Model real RR")
scatter!(uvdist, imag.(getindex.(vis, 1, 2)), color=:yellow, alpha=0.5, marker=:square, label="Model imag RR", legend=false)

p4 = scatter(uvdist, real.(getindex.(vis, 2, 2)), yerr=getindex.(err, 2, 2), title="LL", color=:blue, marker=:circ, label="Data real RR")
scatter!(uvdist, imag.(getindex.(vis, 2, 2)), yerr=getindex.(err, 2, 2), color=:orange, marker=:circ, label="Data imag RR")
scatter!(uvdist, real.(getindex.(vis, 2, 2)), color=:cyan, alpha=0.5, marker=:square, label="Model real RR")
scatter!(uvdist, imag.(getindex.(vis, 2, 2)), color=:yellow, alpha=0.5, marker=:square, label="Model imag RR", legend=false)

plot(p1, p2, p3, p4, layout=(2,2))

# Computing information
# ```
# Julia Version 1.7.3
# Commit 742b9abb4d (2022-05-06 12:58 UTC)
# Platform Info:
#   OS: Linux (x86_64-pc-linux-gnu)
#   CPU: 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-12.0.1 (ORCJIT, tigerlake)
# ```


gAA = CSV.read(joinpath(@__DIR__, "DomTest/gain"))
