function sampling_tau2(iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  obj_prior::PriorsMod1_V6,
  temperature::Float64) where {TD <:GeneralData}
  
  sampling_tau2_cluster(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)


  #if obj_mixture_mcmc.K[1] == obj_mixture_mcmc.Kmax
    
  #  for k in (obj_mixture_mcmc.K[1]+1):obj_mixture_mcmc.Kmax
  #    sampling_tau2_empty(iterations, k, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
  #  end

  #end




end

function sampling_tau2_empty(iterations::Int64,
  k::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::GpData_Vers9,
  obj_data_prop::GpData_Vers9,
  obj_prior::PriorsMod1_V6,
  temperature::Float64)

  obj_data_mcmc.tau2[k] = rand(obj_prior.tau2)
  obj_data_prop.tau2[k] = obj_data_mcmc.tau2[k]

end



function sampling_tau2_cluster(iterations::Int64,
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


  a_p::Vector{Float64} = ones(Float64, obj_mixture_mcmc.K[1]) .* params(obj_prior.tau2)[1]
  b_p::Vector{Float64} = ones(Float64, obj_mixture_mcmc.K[1]) .* params(obj_prior.tau2)[2]
  for iobs in 1:obj_data_mcmc.n_points
    k = obj_mixture_mcmc.cluster[iobs]
    a_p[k] += 0.5
    b_p[k] += 0.5 * ((obj_data_mcmc.obs[iobs] - obj_data_mcmc.gp[iobs, k])^2.0)/ temperature
    
  end

  
  for k in 1:obj_mixture_mcmc.K[1]
    #println([a_p[k], b_p[k]])
#    println([a_p[k], b_p[k], params(obj_prior.tau2)[3], params(obj_prior.tau2)[4]])
    app = rand(InverseGamma(a_p[k], b_p[k]))
    iii = 0
    while (app < params(obj_prior.tau2)[3] || (app > params(obj_prior.tau2)[4])) & (iii <10000)
      app = rand(InverseGamma(a_p[k], b_p[k]))
      iii += 1
    end
    #if k == 1
    #  println("tau2 ", [a_p[k], b_p[k]])
    #end
    
    if iii == 10000
      #println("HERE tau2")  
      app = rand(Normal(obj_data_mcmc.tau2[k], sample(sd_vector, 1)[1]))
      
      MH_ratio = 0.0
      if (app > params(obj_prior.tau2)[3] && (app < params(obj_prior.tau2)[4]))
          
        MH_ratio += logpdf(InverseGamma(params(obj_prior.tau2)[1], params(obj_prior.tau2)[1]), app) - logpdf(InverseGamma(params(obj_prior.tau2)[1], params(obj_prior.tau2)[1]), obj_data_mcmc.tau2[k])
        for iobs in 1:obj_data_mcmc.n_points
          k_app = obj_mixture_mcmc.cluster[iobs]
          if k == k_app
            MH_ratio += logpdf(Normal(obj_data_mcmc.gp[iobs, k], (app * temperature)^0.5), obj_data_mcmc.obs[iobs]) - logpdf(Normal(obj_data_mcmc.gp[iobs, k], (obj_data_prop.tau2[k]*temperature)^0.5), obj_data_mcmc.obs[iobs])
          end
        end

      else
        MH_ratio = log(0.0)
      end

      

      if rand(Uniform(0.0,1.0))< exp(MH_ratio)
        
      else
        app = obj_data_mcmc.tau2[k]
      end

    end

    obj_data_mcmc.tau2[k] = app
    obj_data_prop.tau2[k] = obj_data_mcmc.tau2[k]
  end


end


#function sampling_tau2_cluster(iterations::Int64,
#  obj_graph_mcmc::GraphCluter_Vers6,
#  obj_graph_prop::GraphCluter_Vers6,
#  obj_mixture_mcmc::TestMixture_V5,
#  obj_mixture_prop::TestMixture_V5,
#  obj_data_mcmc::GpDataMarginalized_Vers9,
#  obj_data_prop::GpDataMarginalized_Vers9,
#  obj_prior::PriorsMod1_V6)


#end