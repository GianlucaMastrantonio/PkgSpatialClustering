function sampling_rho_and_gp(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  obj_prior::PriorsMod1_V6,
  temperature::Float64) where {TD <:GeneralData}
  
  sampling_rho_cluster(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)
  sampling_gp_cluster(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)


  #if obj_mixture_mcmc.K[1] == obj_mixture_mcmc.Kmax
    
  #  for k in (obj_mixture_mcmc.K[1]+1):obj_mixture_mcmc.Kmax
  #    sampling_rho_and_gp_empty(iterations, k, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
  #  end

  #end




end

function sampling_rho_and_gp_empty(iterations::Int64,
  k::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::GpData_Vers9,
  obj_data_prop::GpData_Vers9,
  obj_prior::PriorsMod1_V6,
  temperature::Float64)

  obj_data_mcmc.rho[k] = rand(obj_prior.rho)
  obj_data_prop.rho[k] = obj_data_mcmc.rho[k]

  obj_data_mcmc.sigma_mat[k].data .= temperature *obj_data_mcmc.sigma2[k] .* exp.(-obj_data_mcmc.rho[k] .* obj_data_mcmc.distance_mat)
  cholesky_mat = cholesky(obj_data_mcmc.sigma_mat[k])
  obj_data_mcmc.inv_sigma_mat[k].data .= inv(cholesky_mat)
  obj_data_mcmc.log_det[k] = 0.0
  for i = 1:size(cholesky_mat.L.data, 1)
    obj_data_mcmc.log_det[k] += 2.0 * log(cholesky_mat.L.data[i, i])
  end

  obj_data_prop.rho[k] = obj_data_mcmc.rho[k]

  obj_data_prop.sigma_mat[k].data .= obj_data_mcmc.sigma_mat[k].data
  obj_data_prop.inv_sigma_mat[k] .= obj_data_mcmc.inv_sigma_mat[k]
  obj_data_prop.log_det[k] = obj_data_mcmc.log_det[k]

  obj_data_mcmc.gp[:, k] = cholesky_mat.L*rand(Normal(0.0, 1.0), obj_data_mcmc.n_points) .+ obj_data_mcmc.mu[k]
  obj_data_prop.gp[:, k] = obj_data_mcmc.gp[:, k] 

end



#function sampling_rho_tau2_sigma2_empty(iterations::Int64,
#  k::Int64,
#  obj_graph_mcmc::GraphCluter_Vers6,
#  obj_graph_prop::GraphCluter_Vers6,
#  obj_mixture_mcmc::TestMixture_V5,
#  obj_mixture_prop::TestMixture_V5,
#  obj_data_mcmc::GpData_Vers9,
#  obj_data_prop::GpData_Vers9,
#  obj_prior::PriorsMod1_V6)

#  obj_data_mcmc.rho[k] = rand(obj_prior.rho)
#  obj_data_prop.rho[k] = obj_data_mcmc.rho[k]

#  obj_data_mcmc.tau2[k] = rand(obj_prior.tau2)
#  obj_data_prop.tau2[k] = obj_data_mcmc.tau2[k]

#  obj_data_mcmc.sigma2[k] = rand(obj_prior.sigma2)
#  obj_data_prop.sigma2[k] = obj_data_mcmc.sigma2[k]

#  obj_data_mcmc.sigma_mat[k].data .= obj_data_mcmc.sigma2[k] .* exp.(-obj_data_mcmc.rho[k] .* obj_data_mcmc.distance_mat) + obj_data_mcmc.tau2[k] .*I(obj_data_mcmc.n_points)
#  cholesky_mat = cholesky(obj_data_mcmc.sigma_mat[k])
#  obj_data_mcmc.inv_sigma_mat[k].data .= inv(cholesky_mat)
#  obj_data_mcmc.log_det[k] = 0.0
#  for i = 1:size(cholesky_mat.L.data, 1)
#    obj_data_mcmc.log_det[k] += 2.0 * log(cholesky_mat.L.data[i, i])
#  end

  

#  obj_data_prop.sigma_mat[k].data .= obj_data_mcmc.sigma_mat[k].data
#  obj_data_prop.inv_sigma_mat[k] .= obj_data_mcmc.inv_sigma_mat[k]
#  obj_data_prop.log_det[k] = obj_data_mcmc.log_det[k]

  

#end



function sampling_rho_tau2_sigma2_empty(iterations::Int64,
  k::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::GpDataMarginalized_Vers9,
  obj_data_prop::GpDataMarginalized_Vers9,
  obj_prior::PriorsMod1_V6,
  temperature::Float64)

  obj_data_mcmc.rho[k] = rand(obj_prior.rho)
  obj_data_prop.rho[k] = obj_data_mcmc.rho[k]

  obj_data_mcmc.tau2[k] = rand(obj_prior.tau2)
  obj_data_prop.tau2[k] = obj_data_mcmc.tau2[k]

  obj_data_mcmc.sigma2[k] = rand(obj_prior.sigma2)
  obj_data_prop.sigma2[k] = obj_data_mcmc.sigma2[k]

  obj_data_mcmc.sigma_mat[k].data .= temperature * obj_data_mcmc.sigma2[k] .* exp.(-obj_data_mcmc.rho[k] .* obj_data_mcmc.distance_mat) + temperature * obj_data_mcmc.tau2[k] .* I(obj_data_mcmc.n_points)
  
  cholesky_mat = cholesky(obj_data_mcmc.sigma_mat[k])
  obj_data_mcmc.inv_sigma_mat[k].data .= inv(cholesky_mat)
  obj_data_mcmc.log_det[k] = 0.0
  for i = 1:size(cholesky_mat.L.data, 1)
    obj_data_mcmc.log_det[k] += 2.0 * log(cholesky_mat.L.data[i, i])
  end



  obj_data_prop.sigma_mat[k].data .= obj_data_mcmc.sigma_mat[k].data
  obj_data_prop.inv_sigma_mat[k] .= obj_data_mcmc.inv_sigma_mat[k]
  obj_data_prop.log_det[k] = obj_data_mcmc.log_det[k]



end





function sampling_gp_cluster(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  obj_prior::PriorsMod1_V6,
  temperature::Float64) where {TD<:GeneralData}
  
  sigma_post::Symmetric{Float64,Matrix{Float64}} = Symmetric(zeros(Float64, obj_data_mcmc.n_points, obj_data_mcmc.n_points)) 
  mu_post::Vector{Float64} = zeros(Float64, obj_data_mcmc.n_points)
  for k in 1:obj_mixture_mcmc.K[1]
    
  
    sigma_post.data .= obj_data_mcmc.inv_sigma_mat[k]  
    mu_post .= sum(obj_data_mcmc.inv_sigma_mat[k], dims=2) .* obj_data_mcmc.mu[k]
    for i in 1:obj_data_mcmc.n_points
      k_iobs = obj_mixture_mcmc.cluster[i]
      if k == k_iobs
        sigma_post.data[i,i] += 1.0 / (obj_data_mcmc.tau2[k]*temperature)
        mu_post[i] += obj_data_mcmc.obs[i] / (obj_data_mcmc.tau2[k] * temperature)
      end
    end
    
    sigma_post.data .= inv(cholesky(sigma_post))
    mu_post .= sigma_post * mu_post

    obj_data_mcmc.gp[:,k] = rand(MvNormal(mu_post, sigma_post))
    obj_data_prop.gp[:,k] = obj_data_mcmc.gp[:,k]


  end
  

end



function sampling_rho_cluster(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  obj_prior::PriorsMod1_V6,
  temperature::Float64) where {TD<:GeneralData}

  
  sd_vector::Vector{Float64} = [0.0001, 0.001, 0.01, 0.5,1.0]
  MH_ratio::Float64 = 0.0
  
  

  for k in 1:obj_mixture_mcmc.K[1]
    @toggled_assert obj_data_mcmc.rho[k] >= params(obj_prior.rho)[1]
    @toggled_assert obj_data_mcmc.rho[k] <= params(obj_prior.rho)[2]
    MH_ratio = 0.0
    obj_data_prop.rho[k] = rand(Normal(obj_data_mcmc.rho[k], sample(sd_vector, 1)[1]))
  #  obj_data_prop.rho[k] = obj_data_mcmc.rho[k]
    if (obj_data_prop.rho[k]>params(obj_prior.rho)[1]) && (obj_data_prop.rho[k]<params(obj_prior.rho)[2])

      obj_data_prop.sigma_mat[k].data .= temperature * obj_data_prop.sigma2[k] .* exp.(-obj_data_prop.rho[k] .* obj_data_prop.distance_mat)
      cholesky_mat = cholesky(obj_data_prop.sigma_mat[k])
      obj_data_prop.inv_sigma_mat[k].data .= inv(cholesky_mat)
      obj_data_prop.log_det[k] = 0.0
      for i = 1:size(cholesky_mat.L.data, 1)
        obj_data_prop.log_det[k] += 2.0 * log(cholesky_mat.L.data[i, i])
      end

      MH_ratio += -0.5 * obj_data_prop.log_det[k] - 0.5 * transpose(obj_data_prop.gp[:, k] .- obj_data_prop.mu[k]) * obj_data_prop.inv_sigma_mat[k] * (obj_data_prop.gp[:, k] .- obj_data_prop.mu[k])
      MH_ratio -= -0.5 * obj_data_mcmc.log_det[k] - 0.5 * transpose(obj_data_mcmc.gp[:, k] .- obj_data_mcmc.mu[k]) * obj_data_mcmc.inv_sigma_mat[k] * (obj_data_mcmc.gp[:, k] .- obj_data_mcmc.mu[k])

      #println([obj_data_prop.log_det[k], obj_data_mcmc.log_det[k]])
      #@toggled_assert obj_data_prop.rho[k] == obj_data_mcmc.rho[k]
      @toggled_assert obj_data_prop.sigma2[k] == obj_data_mcmc.sigma2[k]
      #println(obj_data_prop.sigma_mat[k].data[1:3,1:3])
      #println(obj_data_mcmc.sigma_mat[k].data[1:3, 1:3])
      #if iterations >100
      #error("test")
      #end
      #println("k=",k," ", exp(MH_ratio))

      if rand(Uniform(0.0, 1.0)) < exp(MH_ratio)
        
        obj_data_mcmc.rho[k] = obj_data_prop.rho[k]

        obj_data_mcmc.sigma_mat[k].data .= obj_data_prop.sigma_mat[k].data
        obj_data_mcmc.inv_sigma_mat[k] .= obj_data_prop.inv_sigma_mat[k]
        obj_data_mcmc.log_det[k] = obj_data_prop.log_det[k]


      else

        obj_data_prop.rho[k] = obj_data_mcmc.rho[k]

        obj_data_prop.sigma_mat[k].data .= obj_data_mcmc.sigma_mat[k].data
        obj_data_prop.inv_sigma_mat[k] .= obj_data_mcmc.inv_sigma_mat[k]
        obj_data_prop.log_det[k] = obj_data_mcmc.log_det[k]
      end
    else


    end
    #obj_data_prop.rho[k] = obj_data_mcmc.rho[k]

    


    #obj_data_mcmc.sigma_mat[k].data .= obj_data_mcmc.sigmaat[k].data ./ obj_data_mcmc.rho[k]
    #obj_data_mcmc.inv_sigma_mat[k].data .= obj_data_mcmc.inv_sigma_mat[k].data .* obj_data_mcmc.rho[k]
    
    #a_p = params(obj_prior.rho)[1] + 0.5 * obj_data_mcmc.n_points
    #b_p = params(obj_prior.rho)[2] + 0.5 * transpose(obj_data_mcmc.gp[:, k]) * obj_data_mcmc.inv_sigma_mat[k] * obj_data_mcmc.gp[:, k]
    
    #obj_data_mcmc.rho[k] = rand(InverseGamma(a_p, b_p))
    #obj_data_prop.rho[k] = obj_data_mcmc.rho[k]

    #obj_data_mcmc.sigma_mat[k].data .= obj_data_mcmc.sigma_mat[k].data .* obj_data_mcmc.rho[k]
    #obj_data_mcmc.inv_sigma_mat[k].data .= obj_data_mcmc.inv_sigma_mat[k].data ./ obj_data_mcmc.rho[k]

    #obj_data_prop.sigma_mat[k].data .= obj_data_mcmc.sigma_mat[k].data
    #obj_data_prop.inv_sigma_mat[k].data .= obj_data_mcmc.inv_sigma_mat[k].data
  end


end




#function sampling_gp_cluster(iterations::Int64,
#  obj_graph_mcmc::GraphCluter_Vers6,
#  obj_graph_prop::GraphCluter_Vers6,
#  obj_mixture_mcmc::TestMixture_V5,
#  obj_mixture_prop::TestMixture_V5,
#  obj_data_mcmc::GpDataMarginalized_Vers9,
#  obj_data_prop::GpDataMarginalized_Vers9,
#  obj_prior::PriorsMod1_V6)

  


#end

#function sampling_rho_cluster(iterations::Int64,
#  obj_graph_mcmc::GraphCluter_Vers6,
#  obj_graph_prop::GraphCluter_Vers6,
#  obj_mixture_mcmc::TestMixture_V5,
#  obj_mixture_prop::TestMixture_V5,
#  obj_data_mcmc::GpDataMarginalized_Vers9,
#  obj_data_prop::GpDataMarginalized_Vers9,
#  obj_prior::PriorsMod1_V6)





#end