function sampling_sigma2(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  obj_prior::PriorsMod1_V6,
  temperature::Float64) where {TD <:GeneralData}
  
  sampling_sigma2_cluster(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)


  #if obj_mixture_mcmc.K[1] == obj_mixture_mcmc.Kmax
    
  #  for k in (obj_mixture_mcmc.K[1]+1):obj_mixture_mcmc.Kmax
  #    sampling_sigma2_empty(iterations, k, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
  #  end

  #end




end

#function sampling_sigma2_empty(iterations::Int64,
#  k::Int64,
#  obj_graph_mcmc::GraphCluter_Vers6,
#  obj_graph_prop::GraphCluter_Vers6,
#  obj_mixture_mcmc::TestMixture_V5,
#  obj_mixture_prop::TestMixture_V5,
#  obj_data_mcmc::GpData_Vers9,
#  obj_data_prop::GpData_Vers9,
#  obj_prior::PriorsMod1_V6,
#  temperature::Float64)

#  obj_data_mcmc.sigma2[k] = rand(obj_prior.sigma2)
#  obj_data_prop.sigma2[k] = obj_data_mcmc.sigma2[k]

#end



function sampling_sigma2_cluster(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  obj_prior::PriorsMod1_V6,
  temperature::Float64) where {TD<:GeneralData}

  a_p::Float64 = 0.0
  b_p::Float64 = 0.0
  save_sigma::Float64 = 0.0

  sd_vector::Vector{Float64} = [0.0001, 0.001, 0.01, 0.5, 1.0]
  MH_ratio::Float64 = 0.0
  #println("sigma2")
  for k in 1:obj_mixture_mcmc.K[1]
    obj_data_mcmc.sigma_mat[k].data .= obj_data_mcmc.sigma_mat[k].data ./ (obj_data_mcmc.sigma2[k] * temperature)
    obj_data_mcmc.inv_sigma_mat[k].data .= obj_data_mcmc.inv_sigma_mat[k].data .* (obj_data_mcmc.sigma2[k] * temperature)

    obj_data_mcmc.log_det[k] = obj_data_mcmc.log_det[k] - obj_data_mcmc.n_points * log((obj_data_mcmc.sigma2[k] * temperature))
    
    a_p = params(obj_prior.sigma2)[1] + 0.5 * obj_data_mcmc.n_points
    b_p = params(obj_prior.sigma2)[2] + 0.5 * transpose(obj_data_mcmc.gp[:, k] .- obj_data_mcmc.mu[k]) * obj_data_mcmc.inv_sigma_mat[k] * (obj_data_mcmc.gp[:, k] .- obj_data_mcmc.mu[k])/temperature
    
    #obj_data_mcmc.sigma2[k] = rand(Truncated(InverseGamma(a_p, b_p), params(obj_prior.sigma2)[3], params(obj_prior.sigma2)[4]))
    app = rand(InverseGamma(a_p, b_p))
    iii = 0
    while (app < params(obj_prior.sigma2)[3] || (app > params(obj_prior.sigma2)[4])) & (iii < 10000)
      app = rand(InverseGamma(a_p, b_p))
      iii += 1
    end
    #if k == 1
    #  println("sigma2 ", [a_p, b_p, temperature, obj_data_mcmc.log_det[k]])
    #end
    
    if iii == 10000
      #println("HERE sigma2")
      app = rand(Normal(obj_data_mcmc.sigma2[k], sample(sd_vector, 1)[1]))

      MH_ratio = 0.0
      if (app > params(obj_prior.sigma2)[3] && (app < params(obj_prior.sigma2)[4]))

        MH_ratio += logpdf(InverseGamma(params(obj_prior.sigma2)[1], params(obj_prior.sigma2)[1]), app) - logpdf(InverseGamma(params(obj_prior.sigma2)[1], params(obj_prior.sigma2)[1]), obj_data_mcmc.sigma2[k])
        MH_ratio += -0.5 * obj_data_mcmc.n_points * log(app * temperature) + 0.5 * obj_data_mcmc.n_points * log((obj_data_mcmc.sigma2[k] * temperature))
        MH_ratio += (-0.5 / (app*temperature) + 0.5 / (obj_data_mcmc.sigma2[k] * temperature)) * transpose(obj_data_mcmc.gp[:, k] .- obj_data_mcmc.mu[k]) * obj_data_mcmc.inv_sigma_mat[k] * (obj_data_mcmc.gp[:, k] .- obj_data_mcmc.mu[k])

      else
        MH_ratio = log(0.0)
      end
      #println(MH_ratio)

      if rand(Uniform(0.0, 1.0)) < exp(MH_ratio)

      else
        app = obj_data_mcmc.sigma2[k]
      end

    end


    obj_data_mcmc.sigma2[k] = app
    #println(obj_data_mcmc.sigma2[k])
    #if (app > params(obj_prior.sigma2)[3] && (app < params(obj_prior.sigma2)[4]))
      
    #end
    #println(obj_data_mcmc.sigma2[k])
    obj_data_prop.sigma2[k] = obj_data_mcmc.sigma2[k]

    
    obj_data_mcmc.sigma_mat[k].data .= obj_data_mcmc.sigma_mat[k].data .* (obj_data_mcmc.sigma2[k] * temperature)

    obj_data_mcmc.inv_sigma_mat[k].data .= obj_data_mcmc.inv_sigma_mat[k].data ./ (obj_data_mcmc.sigma2[k] * temperature)
    
    obj_data_mcmc.log_det[k] = obj_data_mcmc.log_det[k] + obj_data_mcmc.n_points * log(obj_data_mcmc.sigma2[k] * temperature)

    obj_data_prop.sigma_mat[k].data .= obj_data_mcmc.sigma_mat[k].data
    obj_data_prop.inv_sigma_mat[k].data .= obj_data_mcmc.inv_sigma_mat[k].data
    obj_data_prop.log_det[k] = obj_data_mcmc.log_det[k]
  end


end



#function sampling_sigma2_cluster(iterations::Int64,
#  obj_graph_mcmc::GraphCluter_Vers6,
#  obj_graph_prop::GraphCluter_Vers6,
#  obj_mixture_mcmc::TestMixture_V5,
#  obj_mixture_prop::TestMixture_V5,
#  obj_data_mcmc::GpDataMarginalized_Vers9,
#  obj_data_prop::GpDataMarginalized_Vers9,
#  obj_prior::PriorsMod1_V6)



#end