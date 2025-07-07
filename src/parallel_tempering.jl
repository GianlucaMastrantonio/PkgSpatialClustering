

function sampling_mcmc(iterations, obj_graph_mcmc::GraphCluter_Vers6, obj_graph_prop::GraphCluter_Vers6, obj_mixture_mcmc::TestMixture_V5, obj_mixture_prop::TestMixture_V5, obj_data_mcmc::TD, obj_data_prop::TD, obj_prior::PriorsMod1_V6, temperature::Float64) where {TD<:GeneralData}


  #sampling_separator(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, temperature)
  #sampling_separator_jump(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, temperature)
  #sampling_ncluster(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)
  #sampling_w(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)
  println("A")
  ## ! Inizio  LASCIARLO QUI
  @time from_marg_to_full(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)
  @time sampling_mu(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)
  @time sampling_tau2(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)
  @time sampling_sigma2(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)
  @time sampling_rho_and_gp(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)
  # !Fine LASCIARLO QUI
  @time from_full_to_marg(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, temperature)
  
  @time for k in 1:obj_mixture_mcmc.K[1]
    #println("k", k)
    update_param_cluster(obj_data_mcmc, obj_mixture_mcmc, k, temperature)
  end

end



function parallel_tempering(iterations::Int64, obj_graph_mcmc::Vector{GraphCluter_Vers6}, obj_graph_prop::Vector{GraphCluter_Vers6}, obj_mixture_mcmc::Vector{TestMixture_V5}, obj_mixture_prop::Vector{TestMixture_V5}, obj_data_mcmc::Vector{GpDataMarginalized_Vers9}, obj_data_prop::Vector{GpDataMarginalized_Vers9}, obj_prior::PriorsMod1_V6, nsamp::Int64, temperatures::Vector{Float64}, inv_index_temp::Vector{Int64})

  clean_like_p1::Matrix{Float64} = zeros(Float64, size(temperatures, 1), obj_mixture_mcmc[1].Kmax[1])
  clean_like_p2::Matrix{Float64} = zeros(Float64, size(temperatures, 1), obj_mixture_mcmc[1].Kmax[1])

  mh_ratio::Float64 = 0.0


  for it in 1:size(temperatures, 1)
    which_obs = obj_data_mcmc[it].which_obs
    for k in 1:obj_mixture_mcmc[it].K[1]
      
      index = 1:obj_data_mcmc[it].n_which_obs[k]
      w_obs = which_obs[k][index]

      clean_like_p1[it, k] = 0.5 * obj_data_mcmc[it].n_which_obs[k] * log(temperatures[inv_index_temp[it]]) - 0.5 * obj_data_mcmc[it].log_det[k] - 0.5 * obj_data_mcmc[it].n_which_obs[k] * log(2.0 * pi)
      clean_like_p2[it, k] = -0.5 * temperatures[inv_index_temp[it]] * transpose(obj_data_mcmc[it].obs[w_obs] .- obj_data_mcmc[it].mu[k]) * obj_data_mcmc[it].inv_sigma_mat[k][index, index] * (obj_data_mcmc[it].obs[w_obs] .- obj_data_mcmc[it].mu[k])

    end

  end

  

  d1::Int64 = 0
  d2::Int64 = 0
  t1::Int64 = 0
  t2::Int64 = 0
  #s1_index::Int64 = 0
  #s2_index::Int64 = 0
  for isamp in 1:nsamp

    mh_ratio = 0.0
    d1 = rand(1:size(temperatures, 1))
    d2 = rand(1:size(temperatures, 1))
    t1 = inv_index_temp[d1]
    t2 = inv_index_temp[d2]

    if d1 != d2

      for k in 1:obj_mixture_mcmc[d1].K[1]

        mh_ratio += clean_like_p1[d1, k] - 0.5 * obj_data_mcmc[d1].n_which_obs[k] * log(temperatures[t2]) + clean_like_p2[d1, k] / temperatures[t2]
        #mh_ratio += clean_like_p1[s2, k] - 0.5 * obj_data_mcmc[s2].n_which_obs[k] * log(temperatures[s1]) + clean_like_p2[s2][k] / temperatures[s1]

        mh_ratio -= clean_like_p1[d1, k] - 0.5 * obj_data_mcmc[d1].n_which_obs[k] * log(temperatures[t1]) + clean_like_p2[d1, k] / temperatures[t1]
        #mh_ratio -= clean_like_p1[s2, k] - 0.5 * obj_data_mcmc[s2].n_which_obs[k] * log(temperatures[s2]) + clean_like_p2[s2][k] / temperatures[s2]

      end
      for k in 1:obj_mixture_mcmc[d2].K[1]
       
        #mh_ratio += clean_like_p1[s1, k] - 0.5 * obj_data_mcmc[s1].n_which_obs[k] * log(temperatures[s2]) + clean_like_p2[s1][k] / temperatures[s2]
        mh_ratio += clean_like_p1[d2, k] - 0.5 * obj_data_mcmc[d2].n_which_obs[k] * log(temperatures[t1]) + clean_like_p2[d2, k] / temperatures[t1]

        #mh_ratio -= clean_like_p1[s1, k] - 0.5 * obj_data_mcmc[s1].n_which_obs[k] * log(temperatures[s1]) + clean_like_p2[s1][k] / temperatures[s1]
        mh_ratio -= clean_like_p1[d2, k] - 0.5 * obj_data_mcmc[d2].n_which_obs[k] * log(temperatures[t2]) + clean_like_p2[d2, k] / temperatures[t2]

      end
      #println(mh_ratio, " from ", t1, " to ", t2)
      if rand(Uniform(0.0, 1.0)) < exp(mh_ratio)

        inv_index_temp[d1] = t2
        inv_index_temp[d2] = t1

        println("Accept", " from ", t1, " to ", t2)

      end

    end

    

  end

  for it in 1:size(temperatures, 1)

    for k in 1:obj_mixture_mcmc[it].K[1]
      #println("k", k)
      update_param_cluster(obj_data_mcmc[it], obj_mixture_mcmc[it], k, temperatures[inv_index_temp[it]])
      
    end

  end
  
  

end