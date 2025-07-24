
function sampling_all_parameters(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::GpDataMarginalized_Vers9,
  obj_data_prop::GpDataMarginalized_Vers9,
  obj_prior::PriorsMod1_V6,
  temperature::Float64,
  adaptive_mcmc:: Vector{MultiType1_V3})

  
  sd_vector::Vector{Float64} = [0.0001, 0.001, 0.01, 0.5]
  MH_ratio::Float64 = 0.0
  do_samp::Bool = true
  accept::Bool = false
  for k in 1:obj_mixture_mcmc.Kmax

    adaptive_mcmc[k].par_acc[1] =  obj_data_mcmc.sigma2[k]
    adaptive_mcmc[k].par_acc[2] =  obj_data_mcmc.tau2[k]
    adaptive_mcmc[k].par_acc[3] =  obj_data_mcmc.rho[k]
    adaptive_mcmc[k].par_acc[4] =  obj_data_mcmc.mu[k]
    
    sample_amcmc(adaptive_mcmc[k])

    obj_data_prop.sigma2[k] = adaptive_mcmc[k].par_prop[1]
    obj_data_prop.tau2[k] = adaptive_mcmc[k].par_prop[2]
    obj_data_prop.rho[k] = adaptive_mcmc[k].par_prop[3]
    obj_data_prop.mu[k] = adaptive_mcmc[k].par_prop[4]
    
    do_samp = true
    MH_ratio = 0.0
    if (obj_data_prop.sigma2[k] < params(obj_prior.sigma2)[3]) || (obj_data_prop.sigma2[k] > params(obj_prior.sigma2)[4])
      do_samp =  false
      
    end
    
    if (obj_data_prop.rho[k] < params(obj_prior.rho)[1]) || (obj_data_prop.rho[k] > params(obj_prior.rho)[2])
      do_samp =   false
      
    end

    if (obj_data_prop.tau2[k] < params(obj_prior.tau2)[3]) || (obj_data_prop.tau2[k] > params(obj_prior.tau2)[4])
      do_samp =  false
      
    end

    if (obj_data_prop.mu[k] < params(obj_prior.mu)[3]) || (obj_data_prop.mu[k] > params(obj_prior.mu)[4])
      do_samp = false
      
    end

    if do_samp == false
      
      MH_ratio = -Inf

    else

      MH_ratio = logpdf(obj_prior.sigma2, obj_data_prop.sigma2[k]) - logpdf(obj_prior.sigma2, obj_data_mcmc.sigma2[k]) +
                logpdf(obj_prior.rho, obj_data_prop.rho[k]) - logpdf(obj_prior.rho, obj_data_mcmc.rho[k]) +
                logpdf(obj_prior.tau2, obj_data_prop.tau2[k]) - logpdf(obj_prior.tau2, obj_data_mcmc.tau2[k]) +
                logpdf(obj_prior.mu, obj_data_prop.mu[k]) - logpdf(obj_prior.mu, obj_data_mcmc.mu[k])

      if k <= obj_mixture_mcmc.K[1]

        update_param_cluster(obj_data_prop, obj_mixture_mcmc, k, temperature)
        MH_ratio += obj_data_prop.log_likelihood[k] - obj_data_mcmc.log_likelihood[k]

      end


      


    end
    #println(MH_ratio)
    if rand(Uniform(0.0, 1.0)) < exp(MH_ratio)

      copy_from_to(obj_data_prop, obj_data_mcmc, k)
      adaptive_mcmc[k].par_acc .= adaptive_mcmc[k].par_prop
      accept = true
    else

      copy_from_to(obj_data_mcmc, obj_data_prop, k)
      accept = false
    end

    update_amcmc(adaptive_mcmc[k], iterations, accept, min(exp(MH_ratio),1.0))




  end
  



end


