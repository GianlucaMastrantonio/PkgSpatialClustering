function sampling_separator(
  iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  temperature::Float64
) where {TD <:GeneralData}

  
  miss_edge_mcmc = obj_mixture_mcmc.miss_edge
  miss_edge_prop = obj_mixture_prop.miss_edge
  sep_clust_mcmc = obj_mixture_mcmc.sep_clust
  sep_clust_prop = obj_mixture_prop.sep_clust

  n_clust = obj_mixture_mcmc.K[1]
  #n_clust_prop = obj_mixture_prop.K[1]

  do_accept::Bool = false
  z_sample::Int64 = 1
  node_fixed::Int64 = 0
  node_changed::Int64 = 0
  node_proposed::Int64 = 0
  w1_mcmc::Vector{Int64} = [0]
  w2_mcmc::Vector{Int64} = [0]
  w_mcmc::Vector{Int64} = [0]
  w1_prop::Vector{Int64} = [0]
  w2_prop::Vector{Int64} = [0]
  MH_ratio::Float64 = 0.0

  cluster_fixed::Int64 = 0
  cluster_changed::Int64 = 0

  freq_table::Matrix{Int64} = zeros(Float64, n_clust, n_clust)
  best_sum::Float64 = 0.0
  best_perm::Vector{Int64} = zeros(Int64, n_clust)
  for k in 1:(n_clust-1)
    
    # decido quale nodo cambiare
    z_sample = rand(1:2,1)[1]

    node_changed = miss_edge_mcmc[k, z_sample]
    node_fixed = miss_edge_mcmc[k, 3-z_sample]

    #cluster_fixed = sep_clust_mcmc[k, z_sample]
    #cluster_changed = sep_clust_mcmc[k, 3-z_sample]
    # trovo quanti e quali vicini ha nel grafo
    w1_mcmc = obj_graph_mcmc.graph_st[findall(x -> x == node_fixed, obj_graph_mcmc.graph_st[:, 1]), 2]
    w2_mcmc = obj_graph_mcmc.graph_st[findall(x -> x == node_fixed, obj_graph_mcmc.graph_st[:, 2]), 1]
    w_mcmc = unique([w1_mcmc; w2_mcmc])

    # controllo e se non ci siano ripetizioni nei nodi
    #filter!(x -> x != node_changed, w_mcmc)
    for k_app in 1:(n_clust-1)
      if k_app != k
        if miss_edge_mcmc[k_app, 1] == node_fixed
          filter!(x -> x != miss_edge_mcmc[k_app, 2], w_mcmc)
        end
        if miss_edge_mcmc[k_app, 2] == node_fixed
          filter!(x -> x != miss_edge_mcmc[k_app, 1], w_mcmc)
        end
      end
    end

    if size(w_mcmc,1) > 0
      ## proposta
      node_proposed = sample(w_mcmc, 1)[1]
      miss_edge_prop[k, z_sample] = node_proposed
      #sep_clust_prop[k, [1,2]] =  sep_clust_prop[k, [2,1]]
      ## cambio i zeta
      #obj_mixture_prop.cluster[node_fixed] = cluster_changed
      #obj_mixture_prop.cluster[node_proposed] = cluster_fixed
      
      it_worked_zeta = update_zeta(obj_mixture_prop, obj_graph_prop)

      if it_worked_zeta == false
        error("didn't work")
      end
      mapping_mat = matching_cluster!(best_perm, freq_table, obj_mixture_mcmc.cluster, obj_mixture_prop.cluster)
      #println("AA")
      #println(countmap(obj_mixture_mcmc.cluster))
      #println(countmap(obj_mixture_prop.cluster))
      #@assert iterations <10
      map_sep_clust!(mapping_mat, obj_mixture_prop.sep_clust, n_clust)

      for k_app in 1:n_clust
        update_which(obj_data_prop, obj_mixture_prop, k_app)
        update_param_cluster(obj_data_prop, obj_mixture_prop, k_app, temperature)

        update_which(obj_data_mcmc, obj_mixture_mcmc, k_app)
        update_param_cluster(obj_data_mcmc, obj_mixture_mcmc, k_app, temperature)
      end
      

    
    

      # trovo quanti e quali vicini ha nel grafo
      w1_prop = obj_graph_prop.graph_st[findall(x -> x == node_fixed, obj_graph_prop.graph_st[:, 1]), 2]
      w2_prop = obj_graph_prop.graph_st[findall(x -> x == node_fixed, obj_graph_prop.graph_st[:, 2]), 1]
      w_prop = unique([w1_prop; w2_prop])

      # controllo e se non ci siano ripetizioni nei nodi
      #filter!(x -> x != node_proposed, w_prop)
      for k_app in 1:(n_clust-1)
        if k_app != k
          if miss_edge_prop[k_app, 1] == node_fixed
            filter!(x -> x != miss_edge_prop[k_app, 2], w_prop)
          end
          if miss_edge_prop[k_app, 2] == node_fixed
            filter!(x -> x != miss_edge_prop[k_app, 1], w_prop)
          end
        end
        
      end

      MH_ratio = 0.0
      MH_ratio += -log(Float64(size(w_prop, 1))) + log(Float64(size(w_mcmc, 1)))
      for k_app in 1:n_clust
        MH_ratio += obj_data_prop.log_likelihood[k_app]  - obj_data_mcmc.log_likelihood[k_app] 
      end
      
      
      
      #println("MH_ratio=", MH_ratio)
      #println([obj_data_prop.log_likelihood[cluster_changed], obj_data_mcmc.log_likelihood[cluster_changed], obj_data_prop.log_likelihood[cluster_fixed], obj_data_mcmc.log_likelihood[cluster_fixed]])
      if rand(Uniform(0.0,1.0))<exp(MH_ratio)
      
        do_accept = true
      
        
      else
      
        do_accept = false
      end

    end
    #if iterations <= 2
    #  do_accept = true
    #end
    if do_accept
      #  println("accept")
      for k = 1:obj_mixture_prop.Kmax
        copy_from_to(obj_data_prop, obj_data_mcmc, k)
      end
      copy_from_to(obj_mixture_prop, obj_mixture_mcmc)
      #copy_from_to(obj_graph_prop, obj_graph_mcmc)


    else
      #println("reject")
      for k = 1:obj_mixture_prop.Kmax
        copy_from_to(obj_data_mcmc, obj_data_prop, k)
      end
      copy_from_to(obj_mixture_mcmc, obj_mixture_prop)
      #  copy_from_to(obj_graph_mcmc, obj_graph_prop)


    end
    


  end
  
  
  
end




function sampling_separator_jump(
  iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  temperature::Float64
) where {TD <:GeneralData}


  miss_edge_mcmc = obj_mixture_mcmc.miss_edge
  miss_edge_prop = obj_mixture_prop.miss_edge
  sep_clust_mcmc = obj_mixture_mcmc.sep_clust
  sep_clust_prop = obj_mixture_prop.sep_clust

  n_clust = obj_mixture_mcmc.K[1]
  #n_clust_prop = obj_mixture_prop.K[1]

  do_accept::Bool = false
  z_sample::Int64 = 1
  node_fixed::Int64 = 0
  node_changed::Int64 = 0
  node_proposed::Int64 = 0
  w1_mcmc::Vector{Int64} = [0]
  w2_mcmc::Vector{Int64} = [0]
  w_mcmc::Vector{Int64} = [0]
  w1_prop::Vector{Int64} = [0]
  w2_prop::Vector{Int64} = [0]
  MH_ratio::Float64 = 0.0

  cluster_fixed::Int64 = 0
  cluster_chaged::Int64 = 0

  freq_table::Matrix{Int64} = zeros(Float64, n_clust, n_clust)
  best_sum::Float64 = 0.0
  best_perm::Vector{Int64} = zeros(Float64, n_clust)
  
  samp_sep::Vector{Int64} = zeros(Int64,2)
  for iihh in 1:10
    for k in 1:(n_clust-1)

      samp_sep = obj_graph_prop.graph_st[sample(1:size(obj_graph_prop.graph_st, 1)), :]
      do_accept = true
      for k_app in 1:(n_clust-1)

        if (samp_sep[1] == miss_edge_prop[k_app, 1]) & ((samp_sep[2] == miss_edge_prop[k_app, 2]))
          do_accept = false
        end
        if (samp_sep[2] == miss_edge_prop[k_app, 1]) & ((samp_sep[1] == miss_edge_prop[k_app, 2]))
          do_accept = false
        end

      end

      if do_accept == true

        miss_edge_prop[k, :] = samp_sep
        it_worked_zeta = update_zeta(obj_mixture_prop, obj_graph_prop)
        if it_worked_zeta == false
          error("didn't work")
        end
        mapping_mat = matching_cluster!(best_perm, freq_table, obj_mixture_mcmc.cluster, obj_mixture_prop.cluster)
        map_sep_clust!(mapping_mat, obj_mixture_prop.sep_clust, n_clust)

        MH_ratio = 0.0
        for k_app in 1:n_clust
          update_which(obj_data_prop, obj_mixture_prop, k_app)
          update_param_cluster(obj_data_prop, obj_mixture_prop, k_app, temperature)

          update_which(obj_data_mcmc, obj_mixture_mcmc, k_app)
          update_param_cluster(obj_data_mcmc, obj_mixture_mcmc, k_app, temperature)


          MH_ratio += obj_data_prop.log_likelihood[k_app]  - obj_data_mcmc.log_likelihood[k_app] 
        end





        #println("MH_ratio=", MH_ratio)
        #println([obj_data_prop.log_likelihood[cluster_changed], obj_data_mcmc.log_likelihood[cluster_changed], obj_data_prop.log_likelihood[cluster_fixed], obj_data_mcmc.log_likelihood[cluster_fixed]])
        if rand(Uniform(0.0, 1.0)) < exp(MH_ratio)

          do_accept = true


        else

          do_accept = false
        end

      end








      if do_accept
        #  println("accept")
        for k = 1:obj_mixture_prop.Kmax
          copy_from_to(obj_data_prop, obj_data_mcmc, k)
        end
        copy_from_to(obj_mixture_prop, obj_mixture_mcmc)
        #copy_from_to(obj_graph_prop, obj_graph_mcmc)


      else
        #println("reject")
        for k = 1:obj_mixture_prop.Kmax
          copy_from_to(obj_data_mcmc, obj_data_prop, k)
        end
        copy_from_to(obj_mixture_mcmc, obj_mixture_prop)
        #  copy_from_to(obj_graph_mcmc, obj_graph_prop)


      end



    end
  end
  



end


