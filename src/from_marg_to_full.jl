

function from_marg_to_full(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::GpData_Vers9,
  obj_data_prop::GpData_Vers9,
  obj_prior::PriorsMod1_V6,
  temperature::Float64)


end

function from_marg_to_full(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::GpDataMarginalized_Vers9,
  obj_data_prop::GpDataMarginalized_Vers9,
  obj_prior::PriorsMod1_V6,
  temperature::Float64)

   
  ## calcolo le nuove Sigma e simulo
  for k in 1:obj_mixture_mcmc.Kmax
    
    obj_data_mcmc.sigma_mat[k].data .= Symmetric(temperature * obj_data_mcmc.sigma2[k] .* exp.(-obj_data_mcmc.rho[k] .* obj_data_mcmc.distance_mat))

    chol_mat = cholesky(obj_data_mcmc.sigma_mat[k])
    obj_data_mcmc.inv_sigma_mat[k].data .= Symmetric(inv(chol_mat))


    obj_data_mcmc.log_det[k] = 0
    for i = 1:size(chol_mat.L.data, 1)
      obj_data_mcmc.log_det[k] += 2.0 * log(chol_mat.L.data[i, i])
    end


    


    


  end

  sampling_gp_cluster_for_marginal(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)
  for k in 1:obj_mixture_mcmc.Kmax


    index = 1:obj_data_mcmc.n_which_obs[k]
    w_obs = obj_data_mcmc.which_obs[k][index]
    obj_data_mcmc.log_likelihood[k] = (-0.5 * obj_data_mcmc.log_det[k] - 0.5 * obj_data_mcmc.n_which_obs[k] * log(2.0 * pi) - 0.5 * transpose(obj_data_mcmc.obs[w_obs] .- obj_data_mcmc.mu[k]) * obj_data_mcmc.inv_sigma_mat[k][index, index] * (obj_data_mcmc.obs[w_obs] .- obj_data_mcmc.mu[k]))

    
    copy_from_to(obj_data_mcmc, obj_data_prop,k)

    
  end


  #sampling_gp_cluster(iterations::Int64,
  #  obj_graph_mcmc::GraphCluter_Vers6,
  #  obj_graph_prop::GraphCluter_Vers6,
  #  obj_mixture_mcmc::TestMixture_V5,
  #  obj_mixture_prop::TestMixture_V5,
  #  obj_data_mcmc::GpData_Vers9,
  #  obj_data_prop::GpData_Vers9,
  #  obj_prior::PriorsMod1_V6)
end



function sampling_gp_cluster_for_marginal(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::GpDataMarginalized_Vers9,
  obj_data_prop::GpDataMarginalized_Vers9,
  obj_prior::PriorsMod1_V6,
  temperature::Float64)

  sigma_post::Symmetric{Float64,Matrix{Float64}} = Symmetric(zeros(Float64, obj_data_mcmc.n_points, obj_data_mcmc.n_points))
  mu_post::Vector{Float64} = zeros(Float64, obj_data_mcmc.n_points)
  for k in 1:obj_mixture_mcmc.K[1]

    sigma_post.data .= obj_data_mcmc.inv_sigma_mat[k]
    mu_post .= sum(obj_data_mcmc.inv_sigma_mat[k], dims=2) .* obj_data_mcmc.mu[k]

    for i in 1:obj_data_mcmc.n_points
      k_iobs = obj_mixture_mcmc.cluster[i]
      if k == k_iobs
        sigma_post.data[i, i] += 1.0 / (obj_data_mcmc.tau2[k] * temperature)
        mu_post[i] += obj_data_mcmc.obs[i] / (obj_data_mcmc.tau2[k] * temperature)
      end
    end

    sigma_post.data .= inv(cholesky(sigma_post))
    mu_post .= sigma_post * mu_post

    obj_data_mcmc.gp[:, k] = rand(MvNormal(mu_post, sigma_post))
    obj_data_prop.gp[:, k] = obj_data_mcmc.gp[:, k]


  end


end




function from_full_to_marg(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::GpData_Vers9,
  obj_data_prop::GpData_Vers9,
  obj_prior::PriorsMod1_V6,
  temperature::Float64)


end

function from_full_to_marg(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::GpDataMarginalized_Vers9,
  obj_data_prop::GpDataMarginalized_Vers9,
  obj_prior::PriorsMod1_V6,
  temperature::Float64)


  ## calcolo le nuove Sigma e simulo
  for k in 1:obj_mixture_mcmc.Kmax

    obj_data_mcmc.sigma_mat[k].data .= Symmetric(temperature * obj_data_mcmc.sigma2[k] .* exp.(-obj_data_mcmc.rho[k] .* obj_data_mcmc.distance_mat) + temperature *  obj_data_mcmc.tau2[k] .* I(obj_data_mcmc.n_points))

    chol_mat = cholesky(obj_data_mcmc.sigma_mat[k])
    obj_data_mcmc.inv_sigma_mat[k].data .= Symmetric(inv(chol_mat))


    obj_data_mcmc.log_det[k] = 0
    for i = 1:size(chol_mat.L.data, 1)
      obj_data_mcmc.log_det[k] += 2 * log(chol_mat.L.data[i, i])
    end





  end


  for k in 1:obj_mixture_mcmc.Kmax


    copy_from_to(obj_data_mcmc, obj_data_prop, k)


  end


end
