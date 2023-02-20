using Pkg
Pkg.activate(@__DIR__)
using Comrade
using BenchmarkTools
using Distributions
using DistributionsAD
using VLBIImagePriors
using Zygote


using Comrade


# Now we will load some synthetic polarized data.
obs = load_ehtim_uvfits(joinpath(dirname(pathof(Comrade)), "..", "examples", "PolarizedExamples/polarized_gaussian_nogains_withdterms_withfr.uvfits"),
                        joinpath(dirname(pathof(Comrade)), "..", "examples", "PolarizedExamples/array.txt"))
# Notice that, unlike other non-polarized tutorials, we need to include a second argument.
# This is the **array file** of the observation and is required to determine the feed rotation
# of the array.

# Now we scan average the data since the data to boost the SNR and reduce the total data volume.
obs = scan_average(obs)
#-
# Now we extract our observed/corrupted coherency matrices.
dvis = extract_coherency(obs)

function model(θ, metadata)
    (;c, f, p, angparams, dRx, dRy, dLx, dLy, lgp, gpp, lgr, gpr) = θ
    (; grid, cache, tcache, scancache, trackcache) = metadata
    ## Construct the image model
    ## produce Stokes images from parameters
    imgI = f*c
    ## Converts from poincare sphere parameterization of polzarization to Stokes Parameters
    pimg = PoincareSphere2Map(imgI, p, angparams, grid)
    m = ContinuousImage(pimg, cache)

    ## Now construct the basis transformation cache
    jT = jonesT(tcache)

    ## Gain product parameters
    gP = exp.(lgp/2 .+ 1im.*gpp/2)
    Gp = jonesG(gP, gP, scancache)
    ## Gain ratio
    gR = exp.(lgr/2 .+ 1im.*gpr/2)
    Gr = jonesG(gR, inv.(gR), trackcache)
    ##D-terms
    D = jonesD(complex.(dRx, dRy), complex.(dLx, dLy), trackcache)
    ## sandwich all the jones matrices together
    J = Gp*Gr*D*jT
    ## form the complete Jones or RIME model. We use tcache here
    ## to set the reference basis of the model.
    return JonesModel(J, m, tcache)
end



function bench(npix, dvis)
    fovx = μas2rad(50.0)
    fovy = μas2rad(50.0)
    nx = npix
    ny = npix
    grid = imagepixels(fovx, fovy, nx, ny) # image grid
    buffer = IntensityMap(zeros(nx, ny), grid) # buffer to store temporary image
    pulse = BSplinePulse{3}() # pulse we will be using
    cache = create_cache(NFFTAlg(dvis), buffer, pulse) # cache to define the NFFT transform

    tcache = TransformCache(dvis; add_fr=true, ehtim_fr_convention=false)
    #-
    # Next we define our cache that maps quantities e.g., gain products, that change from scan-to-scan.
    scancache = jonescache(dvis, ScanSeg())
    #-
    # Finally, we define our cache that maps quantities, e.g., gain ratios and d-terms, that are constant
    # across a observation night, and we collect everything together.
    trackcache = jonescache(dvis, TrackSeg())
    metadata = (;cache, grid, tcache, scancache, trackcache)

    distamp = (AA = Normal(0.0, 0.1),
           AP = Normal(0.0, 0.1),
           LM = Normal(0.0, 0.1),
           AZ = Normal(0.0, 0.1),
           JC = Normal(0.0, 0.1),
           PV = Normal(0.0, 0.1),
           SM = Normal(0.0, 0.1),
           )
    #-
    # For the phases, we assume that the atmosphere effectively scrambles the gains.
    # Since the gain phases are periodic, we also use broad von Mises priors for all stations.
    distphase = (AA = DiagonalVonMises(0.0, inv(1e-6)),
                 AP = DiagonalVonMises(0.0, inv(π^2)),
                 LM = DiagonalVonMises(0.0, inv(π^2)),
                 AZ = DiagonalVonMises(0.0, inv(π^2)),
                 JC = DiagonalVonMises(0.0, inv(π^2)),
                 PV = DiagonalVonMises(0.0, inv(π^2)),
                 SM = DiagonalVonMises(0.0, inv(π^2)),
               )
    #-
    # However, we can now also use a little additional information about the phase offsets
    # where in most cases, they are much better behaved than the products
    distphase_ratio = (AA = DiagonalVonMises(0.0, inv(1e-6)),
                 AP = DiagonalVonMises(0.0, inv(0.1^2)),
                 LM = DiagonalVonMises(0.0, inv(0.1^2)),
                 AZ = DiagonalVonMises(0.0, inv(0.1^2)),
                 JC = DiagonalVonMises(0.0, inv(0.1^2)),
                 PV = DiagonalVonMises(0.0, inv(0.1^2)),
                 SM = DiagonalVonMises(0.0, inv(0.1^2)),
               )


    # Moving onto the d-terms, here we directly parameterize the real and complex components
    # of the d-terms since they are expected to be complex numbers near the origin. To help enforce
    # this smallness, a weakly informative Normal prior is used.
    distD = ( AA = Normal(0.0, 0.1),
              AP = Normal(0.0, 0.1),
              LM = Normal(0.0, 0.1),
              AZ = Normal(0.0, 0.1),
              JC = Normal(0.0, 0.1),
              PV = Normal(0.0, 0.1),
              SM = Normal(0.0, 0.1),
            )

    prior = (
          c = ImageDirichlet(1.0, nx, ny),
          f = Uniform(0.7, 1.2),
          p = ImageUniform(nx, ny),
          angparams = ImageSphericalUniform(nx, ny),
          dRx = CalPrior(distD, trackcache),
          dRy = CalPrior(distD, trackcache),
          dLx = CalPrior(distD, trackcache),
          dLy = CalPrior(distD, trackcache),
          lgp = CalPrior(distamp, scancache),
          gpp = CalPrior(distphase, scancache),
          lgr = CalPrior(distamp, trackcache),
          gpr = CalPrior(distphase_ratio,trackcache),
          )

    lklhd = RadioLikelihood(model, metadata, dvis)
    post = Posterior(lklhd, prior)

    tpost = asflat(post)

    ℓ = logdensityof(tpost)
    x0 = prior_sample(tpost)

    ftime = @belapsed $(ℓ)($x0)
    gztime = @belapsed Zygote.gradient($ℓ, $x0)
    # gftime = @belapsed AD.gradient($(AD.ForwardDiffBackend{npix}()), $ℓ, $x0)
    return (;npix, ftime=ftime, gradtime = gztime, gradfrac = gztime/ftime)#, gftime, gffime/ftime
end

function bench_noe(npix, dvis)
    fovx = μas2rad(50.0)
    fovy = μas2rad(50.0)
    nx = npix
    ny = npix
    grid = imagepixels(fovx, fovy, nx, ny) # image grid
    buffer = IntensityMap(zeros(nx, ny), grid) # buffer to store temporary image
    pulse = BSplinePulse{3}() # pulse we will be using
    cache = create_cache(NFFTAlg(dvis), buffer, pulse) # cache to define the NFFT transform

    tcache = TransformCache(dvis; add_fr=true, ehtim_fr_convention=false)
    #-
    # Next we define our cache that maps quantities e.g., gain products, that change from scan-to-scan.
    scancache = jonescache(dvis, ScanSeg())
    #-
    # Finally, we define our cache that maps quantities, e.g., gain ratios and d-terms, that are constant
    # across a observation night, and we collect everything together.
    trackcache = jonescache(dvis, TrackSeg())
    metadata = (;cache, grid, tcache, scancache, trackcache)

    distamp = (AA = Normal(0.0, 0.1),
           AP = Normal(0.0, 0.1),
           LM = Normal(0.0, 0.1),
           AZ = Normal(0.0, 0.1),
           JC = Normal(0.0, 0.1),
           PV = Normal(0.0, 0.1),
           SM = Normal(0.0, 0.1),
           )
    #-
    # For the phases, we assume that the atmosphere effectively scrambles the gains.
    # Since the gain phases are periodic, we also use broad von Mises priors for all stations.
    distphase = (AA = DiagonalVonMises(0.0, inv(1e-6)),
                 AP = DiagonalVonMises(0.0, inv(π^2)),
                 LM = DiagonalVonMises(0.0, inv(π^2)),
                 AZ = DiagonalVonMises(0.0, inv(π^2)),
                 JC = DiagonalVonMises(0.0, inv(π^2)),
                 PV = DiagonalVonMises(0.0, inv(π^2)),
                 SM = DiagonalVonMises(0.0, inv(π^2)),
               )
    #-
    # However, we can now also use a little additional information about the phase offsets
    # where in most cases, they are much better behaved than the products
    distphase_ratio = (AA = DiagonalVonMises(0.0, inv(1e-6)),
                 AP = DiagonalVonMises(0.0, inv(0.1^2)),
                 LM = DiagonalVonMises(0.0, inv(0.1^2)),
                 AZ = DiagonalVonMises(0.0, inv(0.1^2)),
                 JC = DiagonalVonMises(0.0, inv(0.1^2)),
                 PV = DiagonalVonMises(0.0, inv(0.1^2)),
                 SM = DiagonalVonMises(0.0, inv(0.1^2)),
               )


    # Moving onto the d-terms, here we directly parameterize the real and complex components
    # of the d-terms since they are expected to be complex numbers near the origin. To help enforce
    # this smallness, a weakly informative Normal prior is used.
    distD = ( AA = Normal(0.0, 0.1),
              AP = Normal(0.0, 0.1),
              LM = Normal(0.0, 0.1),
              AZ = Normal(0.0, 0.1),
              JC = Normal(0.0, 0.1),
              PV = Normal(0.0, 0.1),
              SM = Normal(0.0, 0.1),
            )

    prior = (
          c = reshape(DistributionsAD.TuringDirichlet(1.0, nx*ny), nx, ny),
          f = Uniform(0.7, 1.2),
          p = ImageUniform(nx, ny),
          angparams = ImageSphericalUniform(nx, ny),
          dRx = CalPrior(distD, trackcache),
          dRy = CalPrior(distD, trackcache),
          dLx = CalPrior(distD, trackcache),
          dLy = CalPrior(distD, trackcache),
          lgp = CalPrior(distamp, scancache),
          gpp = CalPrior(distphase, scancache),
          lgr = CalPrior(distamp, trackcache),
          gpr = CalPrior(distphase_ratio,trackcache),
          )

    lklhd = RadioLikelihood(model, metadata, dvis)
    post = Posterior(lklhd, prior)

    tpost = asflat(post)

    ℓ = logdensityof(tpost)
    x0 = prior_sample(tpost)

    ftime = @belapsed $(ℓ)($x0)
    gztime = @belapsed Zygote.gradient($ℓ, $x0)
    # gftime = @belapsed AD.gradient($(AD.ForwardDiffBackend{npix}()), $ℓ, $x0)
    return (;npix, ftime=ftime, gradtime = gztime, gradfrac = gztime/ftime)#, gftime, gffime/ftime
end

using Bijectors
HypercubeTransform.asflat(::DistributionsAD.TuringDirichlet) = SimplexBijector()

import TransformVariables as TV
TV.transform_with

bench(6, dvis)

btimes = [bench(npix, dvis) for npix in 4:2:32]

using DataFrames
using CSV

btimes_struct = DataFrame(btimes)
CSV.write("gradient_benchmark_withenzyme.csv", btimes_struct)

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
