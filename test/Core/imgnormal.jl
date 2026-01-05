using Distributions

@testset "ImgNormalData" begin
    # Load test data to get a working model
    _, vis, amp, lcamp, cphase = load_data()

    # Create a simple image grid
    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 32, 32)
    skym = SkyModel(test_model, test_prior(), g)

    # Create a base posterior to sample from
    base_post = VLBIPosterior(skym, vis)

    # Test with a simple reduction function (total flux)
    @testset "Total Flux Measurement" begin
        # Create reduction function for total flux
        reduction = flux

        # Get a sample image to create mock measurement
        x0 = prior_sample(base_post)
        m0 = skymodel(base_post, x0)
        img0 = intensitymap(m0, g)

        # Create mock measurement with noise
        true_flux = flux(img0)
        measurement = [true_flux]
        noise = [0.1]

        # Create ImgNormalData object
        imgdata = ImgNormalData(reduction, measurement, noise)

        # Test that fields are correctly stored
        @test imgdata.reduction === reduction
        @test imgdata.measurement == measurement
        @test imgdata.noise == noise

        # Create likelihood from the image data
        ℓ = Comrade.makelikelihood(imgdata)

        # Test that we can evaluate the likelihood
        @test ℓ isa Comrade.ConditionedLikelihood

        # Test likelihood evaluation with the image
        lp = logdensityof(ℓ, img0)
        @test lp isa Real
        @test isfinite(lp)

        # Test that likelihood is higher when image flux matches measurement
        # Create image with exact flux
        img_exact = img0 .* (measurement[1] / flux(img0))
        lp_exact = logdensityof(ℓ, img_exact)

        # The exact match should have higher likelihood
        @test lp_exact >= lp
    end

    @testset "Centroid Measurement" begin
        # Create reduction function for centroid
        reduction = img -> begin
            cx, cy = centroid(img)
            return [cx, cy]
        end

        # Get a sample image
        x0 = prior_sample(base_post)
        m0 = skymodel(base_post, x0)
        img0 = intensitymap(m0, g)

        # Create mock measurement
        cx0, cy0 = centroid(img0)
        measurement = [cx0, cy0]
        noise = [μas2rad(5.0), μas2rad(5.0)]

        # Create ImgNormalData object
        imgdata = ImgNormalData(reduction, measurement, noise)

        # Create likelihood
        ℓ = Comrade.makelikelihood(imgdata)

        # Test likelihood evaluation
        lp = logdensityof(ℓ, img0)
        @test lp isa Real
        @test isfinite(lp)

        # Test with a different image
        x1 = prior_sample(base_post)
        m1 = skymodel(base_post, x1)
        img1 = intensitymap(m1, g)
        lp1 = logdensityof(ℓ, img1)
        @test lp1 isa Real
        @test isfinite(lp1)
    end

    @testset "Multiple Measurements" begin
        # Create reduction function that returns multiple quantities
        reduction = img -> begin
            f = flux(img)
            cx, cy = centroid(img)
            return [f, cx, cy]
        end

        # Get a sample image
        x0 = prior_sample(base_post)
        m0 = skymodel(base_post, x0)
        img0 = intensitymap(m0, g)

        # Create mock measurements
        f0 = flux(img0)
        cx0, cy0 = centroid(img0)
        measurement = [f0, cx0, cy0]
        noise = [0.1, μas2rad(5.0), μas2rad(5.0)]

        # Create ImgNormalData object
        imgdata = ImgNormalData(reduction, measurement, noise)

        # Create likelihood
        ℓ = Comrade.makelikelihood(imgdata)

        # Test likelihood evaluation
        lp = logdensityof(ℓ, img0)
        @test lp isa Real
        @test isfinite(lp)

        # Test that the likelihood is composed correctly
        # The reduction should return a vector
        red_result = reduction(img0)
        @test length(red_result) == 3
        @test red_result[1] ≈ f0
        @test red_result[2] ≈ cx0
        @test red_result[3] ≈ cy0
    end

    @testset "Combining ImgNormalData with VLBI data" begin
        # Create image measurement data
        reduction = flux
        measurement = [1.0]
        noise = [0.1]
        imgdata = ImgNormalData(reduction, measurement, noise)

        # Create posterior with both VLBI and image data using keyword argument
        post = VLBIPosterior(skym, vis; imgdata = (imgdata,))

        # Test that we can evaluate the posterior
        x_prior = prior_sample(post)
        lp_post = logdensityof(post, x_prior)
        @test lp_post isa Real
        @test isfinite(lp_post)

        # Test dataproducts
        dp = dataproducts(post)
        @test dp isa Tuple
        @test vis in dp
    end

    @testset "ImgNormal ContinuousImage" begin
        imgdata = ImgNormalData(centroid, [0.0, 0.0], [1e-12, 1e-12])
        g = imagepixels(μas2rad(150.0), μas2rad(150.0), 32, 32)
        skym = SkyModel(testimg, testimg_prior(g), g; metadata = (grid = g,))
        
        post = VLBIPosterior(skym, vis; imgdata = (imgdata,))
        x_prior = prior_sample(post)
        @inferred logdensityof(post, x_prior)
        
        skym = SkyModel(testimg_add, testimg_add_prior(g), g; metadata = (grid = g,))
        post = VLBIPosterior(skym, vis; imgdata = (imgdata,))
        x_prior = prior_sample(post)
        @inferred logdensityof(post, x_prior)


    end
end
