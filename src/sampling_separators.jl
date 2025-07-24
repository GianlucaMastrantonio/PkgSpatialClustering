function sampling_separator(
  iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  temperature::Float64, obj_cohesion_mcmc::TC, obj_cohesion_prop::TC
) where {TD<:GeneralData,TC<:CohesionFunction}

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
    z_sample = rand(1:2, 1)[1]

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

    if size(w_mcmc, 1) > 0
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
        #if iterations == 10
        #  println("k_app ", k_app)
        #  println([obj_data_prop.sigma2[k_app], obj_data_prop.rho[k_app], obj_data_prop.tau2[k_app], obj_data_prop.mu[k_app]])
        #  println(obj_data_prop.sigma_mat[k_app][1:10,1:10])
        #end

        update_which(obj_data_prop, obj_mixture_prop, k_app)
        update_param_cluster(obj_data_prop, obj_mixture_prop, k_app, temperature)


        update_which(obj_data_mcmc, obj_mixture_mcmc, k_app)
        update_param_cluster(obj_data_mcmc, obj_mixture_mcmc, k_app, temperature)
        #println("k_app v2 ", k_app)
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
        MH_ratio += obj_data_prop.log_likelihood[k_app] - obj_data_mcmc.log_likelihood[k_app]
      end
      MH_ratio += compute_cohesion(obj_cohesion_mcmc, obj_mixture_prop) - compute_cohesion(obj_cohesion_mcmc, obj_mixture_mcmc)
      #println("diff= ",compute_cohesion(obj_cohesion_mcmc, obj_mixture_prop) - compute_cohesion(obj_cohesion_mcmc, obj_mixture_mcmc))


      #println("MH_ratio=", MH_ratio)
      #println([obj_data_prop.log_likelihood[cluster_changed], obj_data_mcmc.log_likelihood[cluster_changed], obj_data_prop.log_likelihood[cluster_fixed], obj_data_mcmc.log_likelihood[cluster_fixed]])
      if rand(Uniform(0.0, 1.0)) < exp(MH_ratio)

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
      copy_from_to(obj_graph_prop, obj_graph_mcmc)


    else
      #println("reject")
      for k = 1:obj_mixture_prop.Kmax
        copy_from_to(obj_data_mcmc, obj_data_prop, k)
      end
      copy_from_to(obj_mixture_mcmc, obj_mixture_prop)
      copy_from_to(obj_graph_mcmc, obj_graph_prop)


    end



  end



end




#function sampling_separator_jump(
#  iterations::Int64,
#  obj_graph_mcmc::GraphCluter_Vers6,
#  obj_graph_prop::GraphCluter_Vers6,
#  obj_mixture_mcmc::TestMixture_V5,
#  obj_mixture_prop::TestMixture_V5,
#  obj_data_mcmc::TD,
#  obj_data_prop::TD,
#  temperature::Float64
#) where {TD<:GeneralData}


#  miss_edge_mcmc = obj_mixture_mcmc.miss_edge
#  miss_edge_prop = obj_mixture_prop.miss_edge
#  sep_clust_mcmc = obj_mixture_mcmc.sep_clust
#  sep_clust_prop = obj_mixture_prop.sep_clust

#  n_clust = obj_mixture_mcmc.K[1]
#  #n_clust_prop = obj_mixture_prop.K[1]

#  do_accept::Bool = false
#  z_sample::Int64 = 1
#  node_fixed::Int64 = 0
#  node_changed::Int64 = 0
#  node_proposed::Int64 = 0
#  w1_mcmc::Vector{Int64} = [0]
#  w2_mcmc::Vector{Int64} = [0]
#  w_mcmc::Vector{Int64} = [0]
#  w1_prop::Vector{Int64} = [0]
#  w2_prop::Vector{Int64} = [0]
#  MH_ratio::Float64 = 0.0

#  cluster_fixed::Int64 = 0
#  cluster_chaged::Int64 = 0

#  freq_table::Matrix{Int64} = zeros(Float64, n_clust, n_clust)
#  best_sum::Float64 = 0.0
#  best_perm::Vector{Int64} = zeros(Float64, n_clust)

#  samp_sep::Vector{Int64} = zeros(Int64, 2)
#  for iihh in 1:10
#    for k in 1:(n_clust-1)

#      samp_sep = obj_graph_prop.graph_st[sample(1:size(obj_graph_prop.graph_st, 1)), :]
#      do_accept = true
#      for k_app in 1:(n_clust-1)

#        if (samp_sep[1] == miss_edge_prop[k_app, 1]) & ((samp_sep[2] == miss_edge_prop[k_app, 2]))
#          do_accept = false
#        end
#        if (samp_sep[2] == miss_edge_prop[k_app, 1]) & ((samp_sep[1] == miss_edge_prop[k_app, 2]))
#          do_accept = false
#        end

#      end

#      if do_accept == true

#        miss_edge_prop[k, :] = samp_sep
#        it_worked_zeta = update_zeta(obj_mixture_prop, obj_graph_prop)
#        if it_worked_zeta == false
#          error("didn't work")
#        end
#        mapping_mat = matching_cluster!(best_perm, freq_table, obj_mixture_mcmc.cluster, obj_mixture_prop.cluster)
#        map_sep_clust!(mapping_mat, obj_mixture_prop.sep_clust, n_clust)

#        MH_ratio = 0.0
#        for k_app in 1:n_clust
#          update_which(obj_data_prop, obj_mixture_prop, k_app)
#          update_param_cluster(obj_data_prop, obj_mixture_prop, k_app, temperature)

#          update_which(obj_data_mcmc, obj_mixture_mcmc, k_app)
#          update_param_cluster(obj_data_mcmc, obj_mixture_mcmc, k_app, temperature)


#          MH_ratio += obj_data_prop.log_likelihood[k_app] - obj_data_mcmc.log_likelihood[k_app]
#        end





#        #println("MH_ratio=", MH_ratio)
#        #println([obj_data_prop.log_likelihood[cluster_changed], obj_data_mcmc.log_likelihood[cluster_changed], obj_data_prop.log_likelihood[cluster_fixed], obj_data_mcmc.log_likelihood[cluster_fixed]])
#        if rand(Uniform(0.0, 1.0)) < exp(MH_ratio)

#          do_accept = true


#        else

#          do_accept = false
#        end

#      end








#      if do_accept
#        #  println("accept")
#        for k = 1:obj_mixture_prop.Kmax
#          copy_from_to(obj_data_prop, obj_data_mcmc, k)
#        end
#        copy_from_to(obj_mixture_prop, obj_mixture_mcmc)
#        #copy_from_to(obj_graph_prop, obj_graph_mcmc)


#      else
#        #println("reject")
#        for k = 1:obj_mixture_prop.Kmax
#          copy_from_to(obj_data_mcmc, obj_data_prop, k)
#        end
#        copy_from_to(obj_mixture_mcmc, obj_mixture_prop)
#        #  copy_from_to(obj_graph_mcmc, obj_graph_prop)


#      end



#    end
#  end




#end





function sampling_separator_jump_no_clust_changes(
  iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  temperature::Float64
) where {TD<:GeneralData}



  index_visited::Vector{Int64} = zeros(Int64, obj_graph_mcmc.n_points)
  n_clust = obj_mixture_mcmc.K[1]
  #n_clust_prop = obj_mixture_prop.K[1]


  mh_ratio::Float64 = 0.0

  #Random.seed!(1)
  if n_clust > 1
    obj_graph_prop.weight_mat .= obj_graph_mcmc.weight_mat
    for k in 1:(n_clust-1)
      obj_graph_prop.weight_mat.data[obj_mixture_mcmc.miss_edge[k, 1], obj_mixture_mcmc.miss_edge[k, 2]] = rand(Uniform(0.8, 1.0))
      obj_graph_prop.weight_mat.data[obj_mixture_mcmc.miss_edge[k, 2], obj_mixture_mcmc.miss_edge[k, 1]] = obj_graph_prop.weight_mat.data[obj_mixture_mcmc.miss_edge[k, 1], obj_mixture_mcmc.miss_edge[k, 2]]
    end

    data_sep::Matrix{Int64} = zeros(Int64, sum(1:(n_clust-1)), 2)
    prob_sep::Vector{Float64} = zeros(Float64, size(data_sep, 1))
    n_per_sep::Matrix{Int64} = zeros(Int64, n_clust, n_clust)
    for iobs in 1:size(obj_graph_mcmc.neigh_graph, 1)

      k1 = obj_mixture_mcmc.cluster[iobs]
      for iobs2 in 1:size(obj_graph_mcmc.neigh_graph[iobs], 1)

        k2 = obj_mixture_mcmc.cluster[obj_graph_mcmc.neigh_graph[iobs][iobs2]]

        n_per_sep[k1, k2] += 1
      end

    end

    h::Int64 = 1

    for i1 in 1:(n_clust-1)

      for i2 in (i1+1):n_clust

        if rand(Uniform(0.0, 1.0)) < 0.5
          data_sep[h, 1] = i1
          data_sep[h, 2] = i2
        else
          data_sep[h, 2] = i1
          data_sep[h, 1] = i2
        end

        prob_sep[h] = n_per_sep[i1, i2]

        h += 1

      end

    end


    index_samp::Vector{Int64} = sample(1:size(data_sep, 1), Weights(prob_sep), n_clust - 1, replace=false)

    obj_mixture_prop.sep_clust[1:(n_clust-1), :] = data_sep[index_samp, :]


    k1::Int64 = 0
    k2::Int64 = 0


    is_possible::Bool = true

    is_possible = is_possible & (length(unique(data_sep[index_samp, :][:])) == n_clust)

    for ik in 1:size(index_samp, 1)

      is_possible = is_possible & (n_per_sep[data_sep[index_samp[ik], 1], data_sep[index_samp[ik], 2]] > 0)

    end
    #println([(length(unique(data_sep[index_samp, :][:])) == n_clust), is_possible])
    #println([length(unique(data_sep[index_samp, :][:])), n_clust, is_possible])
    #if iterations == 173
    #  println("is_possible = ", is_possible)
    #end
    if is_possible


      iterator_clust::Int64 = 0
      iterator_row::Int64 = 0
      iterator_col::Int64 = 0
      is_accepted::Bool = true
      sum_iter::Int64 = 0
      for ik in 1:size(index_samp, 1)

        iterator_clust = 0
        iterator_row = 1
        iterator_col = 1
        is_accepted = false
        sum_iter = 0
        #if iterations == 41
        #  #println(countmap(obj_mixture_mcmc.cluster))
        #  println(n_per_sep)
        #  println(data_sep[index_samp, :])
        #end
        while !is_accepted
          sum_iter += 1
          obs1 = iterator_row
          z1 = obj_mixture_mcmc.cluster[obs1]

          #if iterations == 41
          #  println("data sep ", [data_sep[ik, 1], data_sep[ik, 2], sum_iter])
          #end

          #println("a ", [obs1, z1, sum_iter, size(obj_graph_mcmc.neigh_graph[obs1], 1)])
          if z1 == data_sep[index_samp[ik], 1]


            obs2 = obj_graph_mcmc.neigh_graph[obs1][iterator_col]
            z2 = obj_mixture_mcmc.cluster[obs2]
            #println("b ", [obs2, z2, sum_iter])
            #if iterations == 41
            #  println("Z ", [z1,z2, sum_iter])
            #end
            #if iterations == 173
            #  println([ik, iterator_clust, obs1, obs2, obj_mixture_prop.cluster[obs1], obj_mixture_prop.cluster[obs2], obj_mixture_mcmc.cluster[obs1], obj_mixture_mcmc.cluster[obs2]])
            #end
            if z2 == data_sep[index_samp[ik], 2]
              #println("c ", [n_per_sep[z1, z2], iterator_clust, n_per_sep[z1, z2] - iterator_clust])
              if rand(Uniform(0.0, 1.0)) < 1.0 / (n_per_sep[z1, z2] - iterator_clust)

                is_accepted = true

                obj_mixture_prop.miss_edge[ik, 1] = obs1
                obj_mixture_prop.miss_edge[ik, 2] = obs2

              else
                iterator_clust += 1
              end

            end

            if iterator_col == size(obj_graph_mcmc.neigh_graph[obs1], 1)

              iterator_col = 1
              iterator_row += 1

            else

              iterator_col += 1

            end

          else

            iterator_row += 1

          end
          if sum_iter > (size(obj_graph_mcmc.neigh_graph, 1)^2)
            println([iterator_row, iterator_col, iterator_clust])
            error("while non concluse")

          end
        end

      end

      for k in 1:(n_clust-1)
        obj_graph_prop.weight_mat.data[obj_mixture_prop.miss_edge[k, 1], obj_mixture_prop.miss_edge[k, 2]] = rand(Uniform(0.4, 0.6))
        obj_graph_prop.weight_mat.data[obj_mixture_prop.miss_edge[k, 2], obj_mixture_prop.miss_edge[k, 1]] = obj_graph_prop.weight_mat.data[obj_mixture_prop.miss_edge[k, 1], obj_mixture_prop.miss_edge[k, 2]]
      end
      #println("test")
      #println(obj_mixture_prop.miss_edge[1:(n_clust-1), :])
      #println(obj_mixture_mcmc.miss_edge[1:(n_clust-1), :])
      index_visited .= 0
      update_st(obj_graph_prop, index_visited)
      ## controllo se è compatibile con la partizione esistente
      bool_is_the_same = true

      it_worked_zeta = update_zeta(obj_mixture_prop, obj_graph_prop)
      #println("it_worked_zeta = ", it_worked_zeta)
      # controllo che i divisori sono possibili
      for ik = 1:(n_clust-1)

        node_a = obj_mixture_prop.miss_edge[ik, 1]
        node_b = obj_mixture_prop.miss_edge[ik, 2]

        bool_is_the_same = bool_is_the_same & ((node_b in obj_graph_prop.neigh_st[node_a]) || (node_a in obj_graph_prop.neigh_st[node_b]))


      end


      if bool_is_the_same
        if it_worked_zeta == false

          println("zz=", obj_mixture_prop.cluster)
          println("a=", obj_mixture_prop.miss_edge)
          println("b=", obj_mixture_prop.sep_clust)
          error("didn't work 2")
        end
        # controllo zeta
        bool_is_the_same = randindex(obj_mixture_prop.cluster, obj_mixture_mcmc.cluster)[2] == 1

        #println("a = ", obj_mixture_prop.cluster)
        #println("b = ", obj_mixture_mcmc.cluster)
        if bool_is_the_same
          #obj_mixture_mcmc.cluster .= obj_mixture_prop.cluster
          mh_ratio = 0.0
          #MH_ratio += compute_cohesion(obj_cohesion_mcmc, obj_mixture_prop) - compute_cohesion(obj_cohesion_mcmc, obj_mixture_mcmc)
          for ik in 1:(n_clust-1)

            #  mh_ratio += -log(n_per_sep[obj_mixture_mcmc.sep_clust[ik, 1], obj_mixture_mcmc.sep_clust[ik, 1]]) + log(n_per_sep[data_sep[index_samp[ik], 1], data_sep[index_samp[ik], 2]])
            #println(exp(mh_ratio))
          end

          if rand(Uniform(0.0, 1.0)) < exp(mh_ratio)


            # i tre comandi sotto sono giusti, altrimenti si scambiano i cluster 
            obj_mixture_mcmc.miss_edge .= obj_mixture_prop.miss_edge
            obj_mixture_mcmc.sep_clust .= obj_mixture_prop.sep_clust
            obj_mixture_prop.cluster .= obj_mixture_mcmc.cluster
            #copy_from_to(obj_mixture_prop, obj_mixture_mcmc)
            copy_from_to(obj_graph_prop, obj_graph_mcmc)
          else
            copy_from_to(obj_mixture_mcmc, obj_mixture_prop)
            copy_from_to(obj_graph_mcmc, obj_graph_prop)
          end




        else

          #println(obj_mixture_prop.cluster)
          #println(obj_mixture_mcmc.cluster)
          #println(obj_mixture_prop.miss_edge[1:(n_clust-1), :])
          #println(obj_mixture_mcmc.miss_edge[1:(n_clust-1), :])

          #println(n_per_sep)
          #println(data_sep[index_samp, :])


          #@save "/Users/gianlucamastrantonio/Politecnico di Torino Staff Dropbox/Gianluca Mastrantonio/lavori/ClusterSpaziale/real data/Paci_e_Co/analisierrori/SeparatorV2" * string(iterations) * ".jld2" obj_mixture_mcmc obj_mixture_prop n_per_sep data_sep obj_graph_mcmc obj_graph_prop
          println("Didn't work")
          println([length(unique(data_sep[index_samp, :][:])), n_clust])
          println(data_sep[index_samp, :])
          copy_from_to(obj_graph_mcmc, obj_graph_prop)
          copy_from_to(obj_mixture_mcmc, obj_mixture_prop)
          return nothing
          error("Didn't work")

        end

      else

        @save "/Users/gianlucamastrantonio/Politecnico di Torino Staff Dropbox/Gianluca Mastrantonio/lavori/ClusterSpaziale/real data/Paci_e_Co/analisierrori/Separator.jld2" obj_mixture_mcmc obj_mixture_prop n_per_sep data_sep obj_graph_mcmc obj_graph_prop index_samp
        copy_from_to(obj_graph_mcmc, obj_graph_prop)
        copy_from_to(obj_mixture_mcmc, obj_mixture_prop)
        println("Not the same")
        return nothing

        error("Not the same")

      end

    else

      copy_from_to(obj_mixture_mcmc, obj_mixture_prop)
      copy_from_to(obj_graph_mcmc, obj_graph_prop)

    end
    #println(n_per_sep)
    #println(data_sep[index_samp, :])




  end



end


