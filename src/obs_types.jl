abstract type GeneralData end

struct GpData_Vers9 <: GeneralData

    n_points::Int64

    obs::Vector{Float64}
    mu::Vector{Float64}
    sigma2::Vector{Float64}
    tau2::Vector{Float64}
    rho::Vector{Float64}
    gp::Matrix{Float64}
    log_likelihood::Vector{Float64}

    sigma_mat::Vector{Symmetric{Float64,Matrix{Float64}}}
    inv_sigma_mat::Vector{Symmetric{Float64,Matrix{Float64}}}


    coords::Matrix{Float64}
    distance_mat::Symmetric{Float64,Matrix{Float64}}

    which_obs::Vector{Vector{Int64}}
    n_which_obs::Vector{Int64}

    log_det::Vector{Float64}


    function GpData_Vers9(
        n_points::Int64, obs::Vector{Float64}, mu::Vector{Float64}, sigma2::Vector{Float64}, tau2::Vector{Float64}, rho::Vector{Float64}, gp::Matrix{Float64}, log_likelihood::Vector{Float64}, sigma_mat::Vector{Symmetric{Float64,Matrix{Float64}}}, inv_sigma_mat::Vector{Symmetric{Float64,Matrix{Float64}}}, coords::Matrix{Float64}, distance_mat::Symmetric{Float64,Matrix{Float64}}, which_obs::Vector{Vector{Int64}}, n_which_obs::Vector{Int64}, log_det::Vector{Float64}
    )

        new(n_points, obs, mu, sigma2, tau2, rho, gp, log_likelihood, sigma_mat, inv_sigma_mat, coords, distance_mat, which_obs, n_which_obs, log_det)

    end

end



struct GpDataMarginalized_Vers9 <: GeneralData

    n_points::Int64

    obs::Vector{Float64}
    mu::Vector{Float64}
    sigma2::Vector{Float64}
    tau2::Vector{Float64}
    rho::Vector{Float64}
    gp::Matrix{Float64}
    log_likelihood::Vector{Float64}

    sigma_mat::Vector{Symmetric{Float64,Matrix{Float64}}}
    inv_sigma_mat::Vector{Symmetric{Float64,Matrix{Float64}}}


    coords::Matrix{Float64}
    distance_mat::Symmetric{Float64,Matrix{Float64}}

    which_obs::Vector{Vector{Int64}}
    n_which_obs::Vector{Int64}

    log_det::Vector{Float64}


    function GpDataMarginalized_Vers9(
        n_points::Int64, obs::Vector{Float64}, mu::Vector{Float64}, sigma2::Vector{Float64}, tau2::Vector{Float64}, rho::Vector{Float64}, gp::Matrix{Float64}, log_likelihood::Vector{Float64}, sigma_mat::Vector{Symmetric{Float64,Matrix{Float64}}}, inv_sigma_mat::Vector{Symmetric{Float64,Matrix{Float64}}}, coords::Matrix{Float64}, distance_mat::Symmetric{Float64,Matrix{Float64}}, which_obs::Vector{Vector{Int64}}, n_which_obs::Vector{Int64}, log_det::Vector{Float64}
    )

        new(n_points, obs, mu, sigma2, tau2, rho, gp, log_likelihood, sigma_mat, inv_sigma_mat, coords, distance_mat, which_obs, n_which_obs, log_det)

    end

end



function copy_from_to(from::TData, to::TData, k::Int64) where {TData<:GeneralData}
    
    to.mu[k] = from.mu[k]
    to.sigma2[k] = from.sigma2[k]
    to.tau2[k] = from.tau2[k]
    to.rho[k] = from.rho[k]
    #for i in axes(from.gp, 1)
    #    to.gp[i, k] = from.gp[i, k]
    #end
    to.gp[:, k] = from.gp[:, k]

    to.log_likelihood[k] = from.log_likelihood[k]

    #to.sigma_mat[k] = deepcopy(from.sigma_mat[k])
    #to.inv_sigma_mat[k] = deepcopy(from.inv_sigma_mat[k])
    to.sigma_mat[k].data .= from.sigma_mat[k].data
    to.inv_sigma_mat[k].data .= from.inv_sigma_mat[k].data
    
    
    to.which_obs[k] .= from.which_obs[k]
    to.n_which_obs[k] = from.n_which_obs[k]
    to.log_det[k] = from.log_det[k]



end

function GpData_Vers9(; obs::Vector{Float64}, gp_init::Matrix{Float64}, mu_init::Vector{Float64}, sigma2_init::Vector{Float64}, tau2_init::Vector{Float64}, rho_init::Vector{Float64}, obj_mixture::TestMixture_V5, obj_graph::GraphCluter_Vers5)

    n_points = size(obs, 1)

    mu::Vector{Float64} = zeros(Float64, obj_mixture.Kmax[1])
    sigma2::Vector{Float64} = zeros(Float64, obj_mixture.Kmax[1])
    tau2::Vector{Float64} = zeros(Float64, obj_mixture.Kmax[1])
    rho::Vector{Float64} = zeros(Float64, obj_mixture.Kmax[1])
    gp::Matrix{Float64} = zeros(Float64, n_points, obj_mixture.Kmax[1])


    mu .= mu_init
    sigma2 .= sigma2_init
    tau2 .= tau2_init
    rho .= rho_init
    gp .= gp_init

    distance_mat::Symmetric{Float64,Matrix{Float64}} = Symmetric(zeros(Float64, n_points, n_points))
    for i = 1:n_points
        for j = i:n_points
            distance_mat.data[i, j] = sqrt((obj_graph.coords[i, 1] - obj_graph.coords[j, 1])^2 + (obj_graph.coords[i, 2] - obj_graph.coords[j, 2])^2)
        end
    end

    sigma_mat::Vector{Symmetric{Float64,Matrix{Float64}}} = [deepcopy(Symmetric(sigma2[k] .* exp.(-rho[k] .* distance_mat))) for k = 1:obj_mixture.Kmax[1]]




    inv_sigma_mat::Vector{Symmetric{Float64,Matrix{Float64}}} = [deepcopy(Symmetric(inv(cholesky(sigma_mat[k])))) for k = 1:obj_mixture.Kmax[1]]

    log_det::Vector{Float64} = zeros(Float64, obj_mixture.Kmax[1])
    for k in 1:obj_mixture.Kmax[1]
        
        chol_mat = cholesky(sigma_mat[k])
        log_det[k] = 0.0
        for i = 1:size(chol_mat.L.data, 1)
            log_det[k] += 2 * log(chol_mat.L.data[i, i])
        end

    end
    

    log_likelihood::Vector{Float64} = zeros(Float64, obj_mixture.Kmax[1])
    for k = 1:obj_mixture.Kmax[1]
        log_likelihood[k] = -0.5 * n_points * log(2.0 * pi * tau2[k]) - 0.5 * sum((obs .- mu[k]) .^ 2)/ tau2[k]
    end

    which_obs = [zeros(Int64, n_points) for k = 1:obj_mixture.Kmax[1]]
    n_which_obs = zeros(Int64, obj_mixture.Kmax[1])

    out = GpData_Vers9(n_points, obs, mu, sigma2, tau2, rho, gp, log_likelihood, sigma_mat, inv_sigma_mat, obj_graph.coords, distance_mat, which_obs, n_which_obs, log_det)

    for k = 1:obj_mixture.Kmax[1]
        update_which(out, obj_mixture, k)
        update_param_cluster(out, obj_mixture, k)
    end

    return out
    
end




function GpDataMarginalized_Vers9(; obs::Vector{Float64}, gp_init::Matrix{Float64}, mu_init::Vector{Float64}, sigma2_init::Vector{Float64}, tau2_init::Vector{Float64}, rho_init::Vector{Float64}, obj_mixture::TestMixture_V5, obj_graph::GraphCluter_Vers5)::GpDataMarginalized_Vers9

    n_points = size(obs, 1)

    mu::Vector{Float64} = zeros(Float64, obj_mixture.Kmax[1])
    sigma2::Vector{Float64} = zeros(Float64, obj_mixture.Kmax[1])
    tau2::Vector{Float64} = zeros(Float64, obj_mixture.Kmax[1])
    rho::Vector{Float64} = zeros(Float64, obj_mixture.Kmax[1])
    gp::Matrix{Float64} = zeros(Float64, n_points, obj_mixture.Kmax[1])


    mu .= mu_init
    sigma2 .= sigma2_init
    tau2 .= tau2_init
    rho .= rho_init
    gp .= gp_init

    distance_mat::Symmetric{Float64,Matrix{Float64}} = Symmetric(zeros(Float64, n_points, n_points))
    for i = 1:n_points
        for j = i:n_points
            distance_mat.data[i, j] = sqrt((obj_graph.coords[i, 1] - obj_graph.coords[j, 1])^2 + (obj_graph.coords[i, 2] - obj_graph.coords[j, 2])^2)
        end
    end

    sigma_mat::Vector{Symmetric{Float64,Matrix{Float64}}} = [deepcopy(Symmetric(sigma2[k] .* exp.(-rho[k] .* distance_mat))) for k = 1:obj_mixture.Kmax[1]]




    inv_sigma_mat::Vector{Symmetric{Float64,Matrix{Float64}}} = [deepcopy(Symmetric(inv(cholesky(sigma_mat[k])))) for k = 1:obj_mixture.Kmax[1]]

    log_det::Vector{Float64} = zeros(Float64, obj_mixture.Kmax[1])
    for k in 1:obj_mixture.Kmax[1]

        chol_mat = cholesky(sigma_mat[k])
        log_det[k] = 0.0
        for i = 1:size(chol_mat.L.data, 1)
            log_det[k] += 2 * log(chol_mat.L.data[i, i])
        end

    end


    log_likelihood::Vector{Float64} = zeros(Float64, obj_mixture.Kmax[1])
    #for k = 1:obj_mixture.Kmax[1]
    #    log_likelihood[k] = -0.5 * log_det[k] - 0.2*
    #    -0.5 * n_points * log(2.0 * pi * tau2[k]) - 0.5 * sum((obs .- mu[k]) .^ 2) / tau2[k]
    #end

    which_obs = [zeros(Int64, n_points) for k = 1:obj_mixture.Kmax[1]]
    n_which_obs = zeros(Int64, obj_mixture.Kmax[1])

    out = GpDataMarginalized_Vers9(n_points, obs, mu, sigma2, tau2, rho, gp, log_likelihood, sigma_mat, inv_sigma_mat, obj_graph.coords, distance_mat, which_obs, n_which_obs, log_det)

    for k = 1:obj_mixture.K[1]
        update_which(out, obj_mixture, k)
        update_param_cluster(out, obj_mixture, k)
    end

    return out

end

#function compute_log_det_from_chol(k::Int64, cholmat::Vector{Cholesky{Float64,Matrix{Float64}}})::Float64
#    ret::Float64 = 0.0
#    for i = 1:size(cholmat[k].L.data, 1)
#        ret += 2 * log(cholmat[k].L.data[i, i])
#    end
#    return ret
#end

function update_which(obj_data::TData, obj_mixture::TestMixture_V5, k::Int64) where {TData<:GeneralData}

    which_obs = obj_data.which_obs


    obj_data.n_which_obs[k] = 0
    for i = 1:obj_data.n_points
        if obj_mixture.cluster[i] == k

            obj_data.n_which_obs[k] += 1
            which_obs[k][obj_data.n_which_obs[k]] = i

        end
    end

end

function update_param_cluster(obj_data::GpData_Vers9, obj_mixture::TestMixture_V5, k::Int64)
    
    if obj_data.n_which_obs[k] > 0


        #update_which(obj_data::GpData_Vers9, obj_mixture::TestMixture_V5, k::Int64)
        which_obs = obj_data.which_obs


        #obj_data.n_which_obs[k] = 0
        #for i = 1:obj_data.n_points
        #    if obj_mixture.cluster[i] == k

        #        obj_data.n_which_obs[k] += 1
        #        which_obs[k][obj_data.n_which_obs[k]] = i

        #    end
        #end

        #obj_data.sigma_mat[k] = Symmetric(obj_data.sigma2[k] .* exp.(-obj_data.rho[k] .* obj_data.distance_mat[which_obs[k][1:obj_data.n_which_obs[k]], which_obs[k][1:obj_data.n_which_obs[k]]]))
        #obj_data.chol_sigma_mat[k] = cholesky(obj_data.sigma_mat[k])
        #obj_data.inv_sigma_mat[k] = Symmetric(inv(obj_data.chol_sigma_mat[k]))
        #obj_data.log_det[k] = compute_log_det_from_chol(k, obj_data.chol_sigma_mat) 

        obj_data.log_likelihood[k] = 0.0

        for iobs in which_obs[k][1:obj_data.n_which_obs[k]]

            obj_data.log_likelihood[k] += logpdf(Normal( obj_data.gp[iobs, k], obj_data.tau2[k]^0.5), obj_data.obs[iobs])

        end
    else

        obj_data.log_likelihood[k] = log(0.0)

    end
    
    
    
    
    
    
end



function update_param_cluster(obj_data::GpDataMarginalized_Vers9, obj_mixture::TestMixture_V5, k::Int64)
    
    @toggled_assert obj_data.n_which_obs[k] > 0
    #update_which(obj_data::GpData_Vers9, obj_mixture::TestMixture_V5, k::Int64)
    which_obs = obj_data.which_obs


    #obj_data.n_which_obs[k] = 0
    #for i = 1:obj_data.n_points
    #    if obj_mixture.cluster[i] == k

    #        obj_data.n_which_obs[k] += 1
    #        which_obs[k][obj_data.n_which_obs[k]] = i

    #    end
    #end
    index = 1:obj_data.n_which_obs[k]
    w_obs = which_obs[k][index]
    
    
    obj_data.sigma_mat[k].data[index, index] = Symmetric(obj_data.sigma2[k] .* exp.(-obj_data.rho[k] .* obj_data.distance_mat[w_obs, w_obs]) + obj_data.tau2[k] .* I(obj_data.n_which_obs[k]))

    chol_mat = cholesky(obj_data.sigma_mat[k][index, index])
    obj_data.inv_sigma_mat[k].data[index, index] = Symmetric(inv(chol_mat))


    obj_data.log_det[k] = 0
    for i = 1:size(chol_mat.L.data, 1)
        obj_data.log_det[k] += 2 * log(chol_mat.L.data[i, i])
    end

    obj_data.log_likelihood[k] = -0.5 * obj_data.log_det[k] - 0.5 * obj_data.n_which_obs[k] * log(2.0 * pi) - 0.5 * transpose(obj_data.obs[w_obs] .- obj_data.mu[k]) * obj_data.inv_sigma_mat[k][index, index] * (obj_data.obs[w_obs] .- obj_data.mu[k])





    #if obj_data.n_which_obs[k] > 0


        
        
    #else

    #    obj_data.log_likelihood[k] = log(0.0)

    #end






end





function change_label_cluster_and_parameters!(
    obj_data::TD,
    obj_mixture_mcmc::TestMixture_V5,
    cluster::Vector{Int64},
    obj_data_to_copy_from::TD
)::Dict{Int64,Int64} where {TD<:GeneralData}

    Kmax = obj_mixture_mcmc.Kmax[1]
    un_value = sort!(unique(cluster))
    #println(maximum(un_value))
    diff_set = setdiff(1:obj_mixture_mcmc.K[1], un_value)
    #println("diff_set=", diff_set)
    @assert length(diff_set) == 1

    mapping = Dict(zip(un_value, 1:obj_mixture_mcmc.K[1]))
    cluster .= getindex.(Ref(mapping), cluster)

    
    k_sample::Int64 = sample(obj_mixture_mcmc.K[1]:Kmax,1)[1]
    k_sample = obj_mixture_mcmc.K[1]
    ### missing
    

    for k in diff_set[1]:(Kmax-1)

        obj_data.mu[k] = obj_data_to_copy_from.mu[k+1]
        obj_data.sigma2[k] = obj_data_to_copy_from.sigma2[k+1]
        obj_data.tau2[k] = obj_data_to_copy_from.tau2[k+1]
        obj_data.rho[k] = obj_data_to_copy_from.rho[k+1]
        obj_data.gp[:, k] = obj_data_to_copy_from.gp[:, k+1]
        obj_data.log_likelihood[k] = obj_data_to_copy_from.log_likelihood[k+1]
        obj_data.sigma_mat[k].data .= obj_data.sigma_mat[k+1].data
        obj_data.inv_sigma_mat[k].data .= obj_data.inv_sigma_mat[k+1].data

    end


    obj_data.mu[k_sample] = obj_data_to_copy_from.mu[diff_set[1]]
    obj_data.sigma2[k_sample] = obj_data_to_copy_from.sigma2[diff_set[1]]
    obj_data.tau2[k_sample] = obj_data_to_copy_from.tau2[diff_set[1]]
    obj_data.rho[k_sample] = obj_data_to_copy_from.rho[diff_set[1]]
    obj_data.gp[:, k_sample] = obj_data_to_copy_from.gp[:, diff_set[1]]
    obj_data.log_likelihood[k_sample] = obj_data_to_copy_from.log_likelihood[diff_set[1]]
    obj_data.sigma_mat[k_sample].data .= obj_data.sigma_mat[diff_set[1]].data
    obj_data.inv_sigma_mat[k_sample].data .= obj_data.inv_sigma_mat[diff_set[1]].data

    obj_data.mu[Kmax] = obj_data_to_copy_from.mu[k_sample]
    obj_data.sigma2[Kmax] = obj_data_to_copy_from.sigma2[k_sample]
    obj_data.tau2[Kmax] = obj_data_to_copy_from.tau2[k_sample]
    obj_data.rho[Kmax] = obj_data_to_copy_from.rho[k_sample]
    obj_data.gp[:, Kmax] = obj_data_to_copy_from.gp[:, k_sample]
    obj_data.log_likelihood[Kmax] = obj_data_to_copy_from.log_likelihood[k_sample]
    obj_data.sigma_mat[Kmax].data .= obj_data.sigma_mat[k_sample].data
    obj_data.inv_sigma_mat[Kmax].data .= obj_data.inv_sigma_mat[k_sample].data


    return mapping

end


function sampling_from_prior(iterations::Int64,
  k::Int64,
  obj_graph_mcmc::GraphCluter_Vers5,
  obj_graph_prop::GraphCluter_Vers5,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::GpData_Vers9,
  obj_data_prop::GpData_Vers9,
  obj_prior::PriorsMod1_V3)

    
    sampling_rho_and_gp_empty(iterations, k, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
    sampling_tau2_empty(iterations, k, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
    sampling_mu_empty(iterations, k, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
    sampling_sigma2_empty(iterations, k, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)


end



function sampling_from_prior(iterations::Int64,
    k::Int64,
    obj_graph_mcmc::GraphCluter_Vers5,
    obj_graph_prop::GraphCluter_Vers5,
    obj_mixture_mcmc::TestMixture_V5,
    obj_mixture_prop::TestMixture_V5,
    obj_data_mcmc::GpDataMarginalized_Vers9,
    obj_data_prop::GpDataMarginalized_Vers9,
    obj_prior::PriorsMod1_V3)


    sampling_rho_tau2_sigma2_empty(iterations, k, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
    sampling_mu_empty(iterations, k, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
    


end


