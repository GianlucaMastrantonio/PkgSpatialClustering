abstract type GeneralPriors end

struct PriorsMod1_V3 <: GeneralPriors

    
    mu::Truncated{Normal{Float64},Continuous,Float64,Float64,Float64}
    sigma2::Truncated{InverseGamma{Float64},Continuous,Float64,Float64,Float64}
    tau2::Truncated{InverseGamma{Float64},Continuous,Float64,Float64,Float64}
    rho::Uniform{Float64}
    prob::Beta{Float64}
    n_clust::Truncated{Geometric{Float64},Discrete,Float64,Float64,Float64}




    function PriorsMod1_V3(
        mu::Truncated{Normal{Float64},Continuous,Float64,Float64,Float64},
        sigma2::Truncated{InverseGamma{Float64},Continuous,Float64,Float64,Float64},
        tau2::Truncated{InverseGamma{Float64},Continuous,Float64,Float64,Float64},
        rho::Uniform{Float64},
        prob::Beta{Float64},
        n_clust::Truncated{Geometric{Float64},Discrete,Float64,Float64,Float64},
    )

        new(mu, sigma2, tau2, rho, prob, n_clust)

    end

end

function PriorsMod1_V3(;
    mu::Truncated{Normal{Float64},Continuous,Float64,Float64,Float64},
    sigma2::Truncated{InverseGamma{Float64},Continuous,Float64,Float64,Float64},
    tau2::Truncated{InverseGamma{Float64},Continuous,Float64,Float64,Float64},
    rho::Uniform{Float64},
    prob::Beta{Float64},
    n_clust::Geometric{Float64},
    obj_mixture::TestMixture_V5

)

    PriorsMod1_V3(mu, sigma2, tau2, rho, prob, Truncated(n_clust, 0, obj_mixture.Kmax - 1))

end
