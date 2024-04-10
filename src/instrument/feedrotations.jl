function build_feedrotation(obs::EHTArrayConfiguration)

    # read elevation angles for each station
    config = arrayconfig(obs)
    el1 = StructArrays.component(config.data.elevation, 1)
    el2 = StructArrays.component(config.data.elevation, 2)

    # read parallactic angles for each station
    par1 = StructArrays.component(config.data.parallactic, 1)
    par2 = StructArrays.component(config.data.parallactic, 2)


    # get ehtobservation array info
    tarr  = config.tarr
    ants  = tarr.sites
    elevs = tarr.fr_elevation
    pars  = tarr.fr_parallactic
    offs  = tarr.fr_offset


    # get station names
    bls = config.data.baseline
    ant1 = first.(bls)
    ant2 = last.(bls)

    # get multiplicative prefactors
    f_el1  = zero(el1)
    f_par1 = zero(par1)
    f_off1 = zero(el1)
    f_el2  = zero(el2)
    f_par2 = zero(par2)
    f_off2 = zero(el2)
    for i in eachindex(ant1)
        ind1 = findall(==(ant1[i]), ants) |> first
        ind2 = findall(==(ant2[i]), ants) |> first

        f_el1[i]  = elevs[ind1]
        f_el2[i]  = elevs[ind2]

        f_par1[i] = pars[ind1]
        f_par2[i] = pars[ind2]

        f_off1[i] = offs[ind1]
        f_off2[i] = offs[ind2]
    end
    # combine to get field rotations for each station
    FR1 = (f_el1 .* el1) .+ (f_par1 .* par1) .+ f_off1
    FR2 = (f_el2 .* el2) .+ (f_par2 .* par2) .+ f_off2

    if ehtim_fr_convention
        FR1 .*= 2
        FR2 .*= 2
    end
    S = Complex{eltype(FR1)}
    offdiag = fill(zero(eltype(FR1)), length(FR1))
    jF1 = StructArray{SMatrix{2,2,S,4}}((cis.(-FR1), offdiag, offdiag, cis.(FR1)))
    jF2 = StructArray{SMatrix{2,2,S,4}}((cis.(-FR2), offdiag, offdiag, cis.(FR2)))

    return JonesPairs(jF1, jF2)
end
