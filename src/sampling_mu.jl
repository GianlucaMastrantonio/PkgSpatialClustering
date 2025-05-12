function sampling_mu(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers5,
  obj_graph_prop::GraphCluter_Vers5,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  obj_prior::PriorsMod1_V4) where {TD <:GeneralData}
  
  sampling_mu_cluster(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop,  obj_prior)


  #if obj_mixture_mcmc.K[1] == obj_mixture_mcmc.Kmax
    
  #  for k in (obj_mixture_mcmc.K[1]+1):obj_mixture_mcmc.Kmax
  #      sampling_mu_empty(iterations,k,obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
  #  end

  #end




end

function sampling_mu_empty(iterations::Int64,
  k::Int64,
  obj_graph_mcmc::GraphCluter_Vers5,
  obj_graph_prop::GraphCluter_Vers5,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  obj_prior::PriorsMod1_V4) where {TD <:GeneralData}

  obj_data_mcmc.mu[k] = rand(obj_prior.mu)
  obj_data_prop.mu[k] = obj_data_mcmc.mu[k]

end






function sampling_mu_cluster(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers5,
  obj_graph_prop::GraphCluter_Vers5,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  obj_prior::PriorsMod1_V4) where {TD<:GeneralData}

  var_p::Vector{Float64} = ones(Float64, obj_mixture_mcmc.K[1]) .* (1.0 / params(obj_prior.mu)[2]^2.0)
  mean_p::Vector{Float64} = ones(Float64, obj_mixture_mcmc.K[1]) .* (params(obj_prior.mu)[1] / params(obj_prior.mu)[2]^2.0)





  for k in 1:obj_mixture_mcmc.K[1]
    #  @toggled_assert is_pos_def(obj_data_mcmc.inv_sigma_mat[k][1:obj_data_mcmc.n_which_obs[k], 1:obj_data_mcmc.n_which_obs[k]])
    #println(k)
    #println(obj_data_mcmc.which_obs[k])
    #println(obj_data_mcmc.n_which_obs[k])

    #w_obs = obj_data_mcmc.which_obs[k][1:obj_data_mcmc.n_which_obs[k]]
    #var_p[k] += sum(obj_data_mcmc.inv_sigma_mat[k][1:obj_data_mcmc.n_which_obs[k], 1:obj_data_mcmc.n_which_obs[k]])
    #mean_p[k] += sum(transpose(obj_data_mcmc.obs[w_obs]) * obj_data_mcmc.inv_sigma_mat[k][1:obj_data_mcmc.n_which_obs[k], 1:obj_data_mcmc.n_which_obs[k]])

    var_p[k] += sum(obj_data_mcmc.inv_sigma_mat[k])
    mean_p[k] += sum(obj_data_mcmc.inv_sigma_mat[k] * obj_data_mcmc.gp[:, k])

  end

  var_p .= 1.0 ./ var_p
  mean_p .= mean_p .* var_p
  for k in 1:obj_mixture_mcmc.K[1]

    obj_data_mcmc.mu[k] = rand(Truncated(Normal(mean_p[k], sqrt(var_p[k])), params(obj_prior.mu)[3], params(obj_prior.mu)[4]))
    obj_data_prop.mu[k] = obj_data_mcmc.mu[k]

  end


  #var_p::Vector{Float64} = ones(Float64, obj_mixture_mcmc.K[1]) .* (1.0 / params(obj_prior.mu)[2]^2.0)
  #mean_p::Vector{Float64} = ones(Float64, obj_mixture_mcmc.K[1]) .* (params(obj_prior.mu)[1] / params(obj_prior.mu)[2]^2.0)
  #for iobs in 1:obj_data_mcmc.n_points
  #  k = obj_mixture_mcmc.cluster[iobs]
  #  var_p[k] += 1.0/obj_data_mcmc.tau2[k]
  #  mean_p[k] += (obj_data_mcmc.obs[iobs] - obj_data_mcmc.gp[iobs,k]) / obj_data_mcmc.tau2[k]
  #end

  #var_p .= 1.0 ./ var_p
  #mean_p .= mean_p .* var_p
  #for k in 1:obj_mixture_mcmc.K[1]
  #  obj_data_mcmc.mu[k] = rand(Truncated(Normal(mean_p[k], sqrt(var_p[k])), params(obj_prior.mu)[3], params(obj_prior.mu)[4]))
  #  obj_data_prop.mu[k] = obj_data_mcmc.mu[k]
  #end


end

function is_pos_def(M::AbstractMatrix)
  try
    cholesky(M)
    return true
  catch e
    return false
  end
end

#function sampling_mu_cluster(iterations::Int64,
#  obj_graph_mcmc::GraphCluter_Vers5,
#  obj_graph_prop::GraphCluter_Vers5,
#  obj_mixture_mcmc::TestMixture_V5,
#  obj_mixture_prop::TestMixture_V5,
#  obj_data_mcmc::GpDataMarginalized_Vers9,
#  obj_data_prop::GpDataMarginalized_Vers9,
#  obj_prior::PriorsMod1_V4)

#  #var_p::Vector{Float64} = ones(Float64, obj_mixture_mcmc.K[1]) .* (1.0 / params(obj_prior.mu)[2]^2.0)
#  #mean_p::Vector{Float64} = ones(Float64, obj_mixture_mcmc.K[1]) .* (params(obj_prior.mu)[1] / params(obj_prior.mu)[2]^2.0)



  
  
#  #for k in 1:obj_mixture_mcmc.K[1]
#  ##  @toggled_assert is_pos_def(obj_data_mcmc.inv_sigma_mat[k][1:obj_data_mcmc.n_which_obs[k], 1:obj_data_mcmc.n_which_obs[k]])
#  #  #println(k)
#  #  #println(obj_data_mcmc.which_obs[k])
#  #  #println(obj_data_mcmc.n_which_obs[k])
    
#  #  #w_obs = obj_data_mcmc.which_obs[k][1:obj_data_mcmc.n_which_obs[k]]
#  #  #var_p[k] += sum(obj_data_mcmc.inv_sigma_mat[k][1:obj_data_mcmc.n_which_obs[k], 1:obj_data_mcmc.n_which_obs[k]])
#  #  #mean_p[k] += sum(transpose(obj_data_mcmc.obs[w_obs]) * obj_data_mcmc.inv_sigma_mat[k][1:obj_data_mcmc.n_which_obs[k], 1:obj_data_mcmc.n_which_obs[k]])
    
#  #  var_p[k] += sum(obj_data_mcmc.inv_sigma_mat[k])
#  #  mean_p[k] += sum(  obj_data_mcmc.inv_sigma_mat[k] * obj_data_mcmc.gp[:,k])

#  #end

#  #var_p .= 1.0 ./ var_p
#  #mean_p .= mean_p .* var_p
#  #for k in 1:obj_mixture_mcmc.K[1]

#  #  obj_data_mcmc.mu[k] = rand(Truncated(Normal(mean_p[k], sqrt(var_p[k])), params(obj_prior.mu)[3], params(obj_prior.mu)[4]))
#  #  obj_data_prop.mu[k] = obj_data_mcmc.mu[k]

#  #end




#end