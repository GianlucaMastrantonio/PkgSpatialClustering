abstract type CohesionFunction end

struct CohesionFunction1_T1 <: CohesionFunction

    alpha::Vector{Float64}
    
    function CohesionFunction1_T1(alpha::Vector{Float64})

        new(alpha)

    end

end


function CohesionFunction1_T1(;parameter::Float64)

    CohesionFunction1_T1([parameter])

end

function create_same_object(obj::CohesionFunction1_T1)

    return CohesionFunction1_T1(parameter = obj.alpha[1])

end



function compute_cohesion(obj_cohesion::CohesionFunction1_T1, obj_mixture::TestMixture_V5)::Float64

    clust_n = countmap(obj_mixture.cluster)
    #println("clust_n=", clust_n)
    log_ret::Float64 = 0.0
    log_alpha::Float64 = log(obj_cohesion.alpha[1])
    for k in 1:length(clust_n)

        log_ret += log_alpha + 1.0*logfactorial(clust_n[k] -1)

    end
    #println("log_ret=",log_ret)
    return log_ret
end



struct NoCohesionFunction1_T1 <: CohesionFunction



    function NoCohesionFunction1_T1()

        new( )

    end

end

function create_same_object(obj::NoCohesionFunction1_T1)

    return NoCohesionFunction1_T1( )

end

function compute_cohesion(obj_cohesion::NoCohesionFunction1_T1, obj_mixture::TestMixture_V5)::Float64

    
    return 0.0
end