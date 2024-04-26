


function bench(npix, lklhd, dlcamp, )
    fovxy = μas2rad(65.0)
    # Now we can feed in the array information to form the cache. We will be using a DFT since
    # it is efficient for so few pixels
    cache = create_cache(DFTAlg(dlcamp), IntensityMap(rand(npix,npix), fovxy, fovxy, BSplinePulse{3}()))
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


fooE = let lca=dlcamp, cp=dcphase,
           fovx=fovx, fovy=fovy,
           cache = cache,
           U = dlcamp.config.ac[:u],
           V = dlcamp.config.ac[:v],
           plan = cache.plan, phases=cache.phases,
           dmat=dlcamp.config.designmat, dmatc=dcphase.config.designmat
    x->begin
        I = IntensityMap(x, fovx, fovy)
        cI = ContinuousImage(I, BSplinePulse{3}())
        mimg = modelimage(cI, cache)
        vis = visibilitymap(mimg, (U=U, V=V))
        mcp = closure_phases(vis, dcphase.config)
        mlca = logclosure_amplitudes(vis, dlcamp.config)
        l1 = sum(abs2, (lca[:amp] .- mlca)./lca[:error])
        l2 = sum(abs2, (cp[:phase] .- mcp)./cp[:error])
        return -0.5*(l1 + l2)
    end
end
