function compute_proposal_change_cluster(;from_k::Int64, to_k::Int64,  Kmax::Int64, nobs::Int64)::Float64

  prob_proposal::Float64 = 0.0 
  prob_which::Float64 = 0.0
  prob_order::Float64 = 0.0
  prob_node::Float64 = 0.0
  
  if from_k == 1
    prob_proposal = -log(2.0)
    
    

  elseif from_k == Kmax
    prob_proposal = -log(2.0)
  
  else

    prob_proposal = -log(3.0)
  

  end

  return prob_proposal - logfactorial(to_k)

end


function OLD_compute_proposal_change_cluster(; from_k::Int64, to_k::Int64, Kmax::Int64, nobs::Int64)::Float64

  prob_proposal::Float64 = 0.0
  prob_which::Float64 = 0.0
  prob_order::Float64 = 0.0
  prob_node::Float64 = 0.0

  if from_k == 1
    prob_proposal = -log(2.0)

    if from_k == to_k
      prob_which = 0.0
      prob_order = 0.0
      prob_node = 0.0
    end

    if from_k < to_k
      prob_which = -log(Kmax - 1)
      prob_order = -logfactorial(1 + 1)
      prob_node = -log(nobs - from_k)
    end

  elseif from_k == Kmax

    if from_k == to_k
      prob_which = -log(Kmax)
      prob_order = -logfactorial(from_k)
      prob_node = -log(nobs - (from_k - 1))
    end

    if from_k > to_k
      prob_which = -log(from_k - 1)
      prob_order = -logfactorial(from_k - 1)
      #prob_node = -log(nobs - (from_k - 2))
      prob_node = 0.0
    end
  else


    prob_proposal = -log(3.0)
    if from_k > to_k
      prob_which = -log(from_k - 1)
      prob_order = -logfactorial(from_k - 1)
      prob_node = 0.0
    end

    if from_k == to_k
      prob_which = -log(from_k)
      prob_order = -logfactorial(from_k)
      prob_node = -log(nobs - (from_k - 1))
    end

    if from_k < to_k
      prob_which = -log(Kmax - from_k)
      prob_order = -logfactorial(from_k + 1)
      prob_node = -log(nobs - (from_k))

    end

  end

  return prob_proposal + prob_which + prob_order

end

function sampling_ncluster(
  iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers6,
  obj_graph_prop::GraphCluter_Vers6,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD,
  obj_prior::PriorsMod1_V6,
  temperature::Float64
) where {TD<:GeneralData}


  #k_mcmc = obj_mixture_mcmc.K[1]
  Kmax = obj_mixture_mcmc.Kmax

  miss_edge_mcmc = obj_mixture_mcmc.miss_edge
  miss_edge_prop = obj_mixture_prop.miss_edge
  #sep_clust_mcmc = obj_mixture_mcmc.sep_clust
  #sep_clust_prop = obj_mixture_prop.sep_clust

  MH_ratio::Float64 = 0.0


  sd_vector::Vector{Float64} = [0.0001, 0.001, 0.01, 0.5, 1.0]

  do_accept::Bool = false
  n_clust::Int64 = obj_mixture_mcmc.K[1]
  n_clust_prop::Int64 = n_clust
  #log_prob_prop::Float64 = 0.0
  #log_prob_mcmc::Float64 = 0.0


  #prob_proposta_prop::Float64 = 0.0
  #prob_proposta_mcmc::Float64 = 0.0
  #prob_which_prop::Float64 = 0.0
  #prob_which_mcmc::Float64 = 0.0
  #prob_order_prop::Float64 = 0.0
  #prob_order_mcmc::Float64 = 0.0

  is_new::Bool = false
  best_perm::Vector{Int64} = zeros(Float64, Kmax)
  freq_table::Matrix{Int64} = zeros(Float64, Kmax, Kmax)
  for isim in 1:100


    obj_mixture_prop.prob[1] = rand(Normal(obj_mixture_mcmc.prob[1], sample(sd_vector, 1)[1]))

    MH_ratio = 0.0
    if (obj_mixture_prop.prob[1] > 0.0) & (obj_mixture_prop.prob[1] < 1.0)

      MH_ratio += (n_clust + params(obj_prior.prob)[1] - 2.0) * log(obj_mixture_prop.prob[1]) + params(obj_prior.prob)[2] * log(1.0 - obj_mixture_prop.prob[1]) - log(1.0 - obj_mixture_prop.prob[1]^Kmax)

      MH_ratio -= (n_clust + params(obj_prior.prob)[1] - 2.0) * log(obj_mixture_mcmc.prob[1]) + params(obj_prior.prob)[2] * log(1.0 - obj_mixture_mcmc.prob[1]) - log(1.0 - obj_mixture_mcmc.prob[1]^Kmax)

      if rand(Uniform(0.0, 1.0)) < exp(MH_ratio)
        obj_mixture_mcmc.prob[1] = obj_mixture_prop.prob[1]
      else

      end

    end
    obj_mixture_prop.prob[1] = obj_mixture_mcmc.prob[1]
  end
  #println("prob", obj_mixture_mcmc.prob[1])
  for isim in 1:20


    MH_ratio = 0.0

    if n_clust == 1
      n_clust_prop = sample(1:2, 1)[1]
    elseif n_clust == Kmax
      n_clust_prop = sample([Kmax, Kmax - 1], 1)[1]
    else
      n_clust_prop = sample([n_clust - 1, n_clust, n_clust + 1], 1)[1]
    end

    #if n_clust_prop == 1
    #  prob_proposta_mcmc = -log(2.0)
    #  prob_order_mcmc = -logfactorial(n_clust)
    #  prob_which_mcmc = -log(Kmax - n_clust_prop)
    #elseif n_clust_prop == Kmax
    #  prob_proposta_mcmc = -log(2.0)
    #  prob_order_mcmc = -logfactorial(n_clust)
    #  prob_which_mcmc = -log(n_clust_prop - 1)
    #else
    #  prob_proposta_mcmc = -log(3.0)
    #  prob_order_mcmc = -logfactorial(n_clust)
    #  if n_clust > n_clust_prop
    #    prob_which_mcmc = -log(Kmax - n_clust_prop)
    #  else
    #    prob_which_mcmc = -log(n_clust_prop - 1)
    #  end
    #end



    obj_mixture_prop.K[1] = n_clust_prop

    if n_clust_prop != n_clust
      MH_ratio = 0.0
      obj_mixture_prop.K[1] = n_clust_prop

      if n_clust_prop > n_clust

        is_new = false

        while is_new == false
          #println("is_new", is_new)
          miss_edge_prop[n_clust_prop-1, :] = obj_graph_prop.graph_st[sample(1:size(obj_graph_prop.graph_st, 1))[1], :]
          #println(miss_edge_prop[n_clust_prop-1, :])
          #println(miss_edge_prop[1:(n_clust_prop-2), :])
          is_new = true
          if n_clust_prop > 2
            for ik in 1:(n_clust_prop-2)

              if (miss_edge_prop[n_clust_prop-1, 1] == miss_edge_prop[ik, 1]) & (miss_edge_prop[n_clust_prop-1, 2] == miss_edge_prop[ik, 2])
                is_new = false
              end
              if (miss_edge_prop[n_clust_prop-1, 2] == miss_edge_prop[ik, 1]) & (miss_edge_prop[n_clust_prop-1, 1] == miss_edge_prop[ik, 2])
                is_new = false
              end


            end
          end

        end

        it_worked_zeta = update_zeta(obj_mixture_prop, obj_graph_prop)

        sampling_from_prior(iterations,
          n_clust_prop,
          obj_graph_mcmc,
          obj_graph_prop,
          obj_mixture_mcmc,
          obj_mixture_prop,
          obj_data_mcmc,
          obj_data_prop,
          obj_prior,
          temperature)
        #  println("A")
        #println([sort!(unique(obj_mixture_mcmc.cluster)), sort!(unique(obj_mixture_prop.cluster))])
        mapping_mat = matching_cluster_different_sizes!(best_perm, freq_table, obj_mixture_mcmc.cluster, obj_mixture_prop.cluster)
        map_sep_clust!(mapping_mat, obj_mixture_prop.sep_clust, n_clust)
        #  println([sort!(unique(obj_mixture_mcmc.cluster)), sort!(unique(obj_mixture_prop.cluster))])
        for k in 1:n_clust_prop
          update_which(obj_data_prop, obj_mixture_prop, k)
          update_param_cluster(obj_data_prop, obj_mixture_prop, k, temperature)
        end

        for k in 1:n_clust
          MH_ratio += obj_data_prop.log_likelihood[k]  - obj_data_mcmc.log_likelihood[k]  
        end
        MH_ratio += obj_data_prop.log_likelihood[n_clust_prop] 

        #  println([MH_ratio, (n_clust_prop - 1) * log(obj_mixture_mcmc.prob[1]), (n_clust - 1) * log(obj_mixture_mcmc.prob[1])])
        # priors
        #println("Piu", [MH_ratio, (n_clust_prop - 1) * log(obj_mixture_mcmc.prob[1]) - (n_clust - 1) * log(obj_mixture_mcmc.prob[1]), -compute_proposal_change_cluster(from_k=n_clust, to_k=n_clust_prop, Kmax=Kmax, nobs=obj_data_mcmc.n_points) + compute_proposal_change_cluster(from_k=n_clust_prop, to_k=n_clust, Kmax=Kmax, nobs=obj_data_mcmc.n_points)])

        MH_ratio += (n_clust_prop - 1) * log(obj_mixture_mcmc.prob[1]) - (n_clust - 1) * log(obj_mixture_mcmc.prob[1])
      
        # proposal
        MH_ratio += -compute_proposal_change_cluster(from_k=n_clust, to_k=n_clust_prop, Kmax=Kmax, nobs=obj_data_mcmc.n_points) + compute_proposal_change_cluster(from_k=n_clust_prop, to_k=n_clust, Kmax=Kmax, nobs=obj_data_mcmc.n_points)
        #println("A=", -compute_proposal_change_cluster(from_k=n_clust, to_k=n_clust_prop, Kmax = Kmax) +compute_proposal_change_cluster(from_k=n_clust_prop, to_k=n_clust, Kmax = Kmax))
        



        if rand(Uniform(0.0, 1.0)) < exp(MH_ratio)
          do_accept = true
        else
          do_accept = false
        end


      end

      if n_clust_prop < n_clust
        MH_ratio = 0.0
        if n_clust_prop > 1
          select_remove = sample(1:(n_clust_prop-1), 1)[1]
          #println("select_remove", select_remove)
          if n_clust_prop > 1
            for imiss in select_remove:(n_clust_prop-1)
              miss_edge_prop[imiss, :] = miss_edge_mcmc[imiss+1, :]
            end
          end
        end

        it_worked_zeta = update_zeta(obj_mixture_prop, obj_graph_prop)
        #println("A")
        #println([sort!(unique(obj_mixture_mcmc.cluster)), sort!(unique(obj_mixture_prop.cluster))])
        mapping_mat = matching_cluster_different_sizes!(best_perm, freq_table, obj_mixture_mcmc.cluster, obj_mixture_prop.cluster)
        map_sep_clust!(mapping_mat, obj_mixture_prop.sep_clust, n_clust_prop)


        mapping_mat = change_label_cluster_and_parameters!(obj_data_prop, obj_mixture_mcmc, obj_mixture_prop.cluster, obj_data_mcmc)
        map_sep_clust!(mapping_mat, obj_mixture_prop.sep_clust, n_clust_prop)



        #    println([sort!(unique(obj_mixture_mcmc.cluster)), sort!(unique(obj_mixture_prop.cluster))])

        for k in 1:n_clust_prop
          update_which(obj_data_prop, obj_mixture_prop, k)
          update_param_cluster(obj_data_prop, obj_mixture_prop, k, temperature)
        end

        for k in 1:(n_clust-1)
          MH_ratio += obj_data_prop.log_likelihood[k] - obj_data_mcmc.log_likelihood[k]
        end
        MH_ratio += -obj_data_mcmc.log_likelihood[n_clust]
        #println("Meno", [MH_ratio, (n_clust_prop - 1) * log(obj_mixture_mcmc.prob[1]) - (n_clust - 1) * log(obj_mixture_mcmc.prob[1]), -compute_proposal_change_cluster(from_k=n_clust, to_k=n_clust_prop, Kmax=Kmax, nobs=obj_data_mcmc.n_points) + compute_proposal_change_cluster(from_k=n_clust_prop, to_k=n_clust, Kmax=Kmax, nobs=obj_data_mcmc.n_points)])
        
        # priors


        MH_ratio += (n_clust_prop - 1) * log(obj_mixture_mcmc.prob[1]) - (n_clust - 1) * log(obj_mixture_mcmc.prob[1])

        # proposal
        MH_ratio += -compute_proposal_change_cluster(from_k=n_clust, to_k=n_clust_prop, Kmax=Kmax, nobs=obj_data_mcmc.n_points) + compute_proposal_change_cluster(from_k=n_clust_prop, to_k=n_clust, Kmax=Kmax, nobs=obj_data_mcmc.n_points)
        #println("A=", -compute_proposal_change_cluster(from_k=n_clust, to_k=n_clust_prop, Kmax=Kmax) + compute_proposal_change_cluster(from_k=n_clust_prop, to_k=n_clust, Kmax=Kmax))
        if rand(Uniform(0.0, 1.0)) < exp(MH_ratio)
          do_accept = true
        else
          do_accept = false
        end


      end

    else

      ### sono uguali
      if n_clust_prop >1
        is_new = false
        select_remove = sample(1:(n_clust-1), 1)[1]
        while is_new == false
          #println("is_new", is_new)
          miss_edge_prop[select_remove, :] = obj_graph_prop.graph_st[sample(1:size(obj_graph_prop.graph_st, 1))[1], :]
          #println(miss_edge_prop[n_clust_prop-1, :])
          #println(miss_edge_prop[1:(n_clust_prop-2), :])
          is_new = true
          if n_clust_prop > 2
            for ik in 1:(n_clust-1)
              if select_remove != ik
                if (miss_edge_prop[select_remove, 1] == miss_edge_prop[ik, 1]) & (miss_edge_prop[select_remove, 2] == miss_edge_prop[ik, 2])
                  is_new = false
                end
                if (miss_edge_prop[select_remove, 2] == miss_edge_prop[ik, 1]) & (miss_edge_prop[select_remove, 1] == miss_edge_prop[ik, 2])
                  is_new = false
                end
              end
              


            end
          end

        end

        it_worked_zeta = update_zeta(obj_mixture_prop, obj_graph_prop)

        #sampling_from_prior(iterations,
        #  n_clust_prop,
        #  obj_graph_mcmc,
        #  obj_graph_prop,
        #  obj_mixture_mcmc,
        #  obj_mixture_prop,
        #  obj_data_mcmc,
        #  obj_data_prop,
        #  obj_prior)
        #  println("A")
        #println([sort!(unique(obj_mixture_mcmc.cluster)), sort!(unique(obj_mixture_prop.cluster))])
        mapping_mat = matching_cluster!(best_perm, freq_table, obj_mixture_mcmc.cluster, obj_mixture_prop.cluster)
        map_sep_clust!(mapping_mat, obj_mixture_prop.sep_clust, n_clust)
        #  println([sort!(unique(obj_mixture_mcmc.cluster)), sort!(unique(obj_mixture_prop.cluster))])
        for k in 1:n_clust
          update_which(obj_data_prop, obj_mixture_prop, k)
          update_param_cluster(obj_data_prop, obj_mixture_prop, k, temperature)
        end
        MH_ratio = 0.0
        for k in 1:n_clust
          MH_ratio += obj_data_prop.log_likelihood[k]   - obj_data_mcmc.log_likelihood[k]  
        end
        #  println("Uguale", [MH_ratio])

        #
        #
        #  println([MH_ratio, (n_clust_prop - 1) * log(obj_mixture_mcmc.prob[1]), (n_clust - 1) * log(obj_mixture_mcmc.prob[1])])
        # priors
        #MH_ratio += (n_clust_prop - 1) * log(obj_mixture_mcmc.prob[1]) - (n_clust - 1) * log(obj_mixture_mcmc.prob[1])

        # proposal
        #MH_ratio += -compute_proposal_change_cluster(from_k=n_clust, to_k=n_clust_prop, Kmax=Kmax, nobs=obj_data_mcmc.n_points) + compute_proposal_change_cluster(from_k=n_clust_prop, to_k=n_clust, Kmax=Kmax, nobs=obj_data_mcmc.n_points)
        #println("A=", -compute_proposal_change_cluster(from_k=n_clust, to_k=n_clust_prop, Kmax = Kmax) +compute_proposal_change_cluster(from_k=n_clust_prop, to_k=n_clust, Kmax = Kmax))


        if rand(Uniform(0.0, 1.0)) < exp(MH_ratio)
          do_accept = true
        else
          do_accept = false
        end
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
      copy_from_to(obj_graph_prop, obj_graph_mcmc)

    else

      #println("reject")
      for k = 1:obj_mixture_prop.Kmax
        copy_from_to(obj_data_mcmc, obj_data_prop, k)
      end
      copy_from_to(obj_mixture_mcmc, obj_mixture_prop)
      copy_from_to(obj_graph_mcmc, obj_graph_prop)

    end
    #println("From ", n_clust, " ", n_clust_prop)
    #println(do_accept)
    #println((obj_mixture_mcmc.K[1]))
    #println((miss_edge_mcmc))
    if (obj_mixture_mcmc.K[1]) <= size(miss_edge_mcmc, 1)
      for imiss in (obj_mixture_mcmc.K[1]):size(miss_edge_mcmc, 1)
        miss_edge_mcmc[imiss, :] .= 0
      end
    end
    n_clust = obj_mixture_mcmc.K[1]
    #println("A")
    #println(obj_mixture_mcmc.K[1], " ", obj_mixture_prop.K[1])
    #println(obj_data_mcmc.n_which_obs[1:obj_mixture_mcmc.K[1]], " ", obj_data_prop.n_which_obs[1:obj_mixture_prop.K[1]])
    #println("B")
    #println("Post")
    #println((obj_mixture_mcmc.K[1]))
    #println((miss_edge_mcmc))
  end

end


#function sampling_ncluster(
#  iterations::Int64,
#  obj_graph_mcmc::GraphCluter_Vers6,
#  obj_graph_prop::GraphCluter_Vers6,
#  obj_mixture_mcmc::TestMixture_V5,
#  obj_mixture_prop::TestMixture_V5,
#  obj_data_mcmc::TD,
#  obj_data_prop::TD,
#  obj_prior::PriorsMod1_V6
#) where {TD <:GeneralData}


#  #k_mcmc = obj_mixture_mcmc.K[1]
#  Kmax = obj_mixture_mcmc.Kmax

#  miss_edge_mcmc = obj_mixture_mcmc.miss_edge
#  miss_edge_prop = obj_mixture_prop.miss_edge
#  #sep_clust_mcmc = obj_mixture_mcmc.sep_clust
#  #sep_clust_prop = obj_mixture_prop.sep_clust

#  MH_ratio::Float64 = 0.0
  

#  sd_vector::Vector{Float64} = [0.0001, 0.001, 0.01, 0.5, 1.0]
  
#  do_accept::Bool = false
#  n_clust::Int64 = obj_mixture_mcmc.K[1]
#  n_clust_prop::Int64 = n_clust
#  log_prob_prop::Float64 = 0.0
#  log_prob_mcmc::Float64 = 0.0
#  is_new::Bool = false
#  best_perm::Vector{Int64} = zeros(Float64, Kmax)
#  freq_table::Matrix{Int64} = zeros(Float64, Kmax, Kmax)
#  for isim in 1:100

    
#    obj_mixture_prop.prob[1] = rand(Normal(obj_mixture_mcmc.prob[1], sample(sd_vector, 1)[1]))
    
#    MH_ratio = 0.0
#    if (obj_mixture_prop.prob[1]  >0.0) & (obj_mixture_prop.prob[1] <1.0)

#      MH_ratio += (n_clust + params(obj_prior.prob)[1] - 2.0) * log(obj_mixture_prop.prob[1]) + params(obj_prior.prob)[2] * log(1.0 - obj_mixture_prop.prob[1]) - log(1.0 - obj_mixture_prop.prob[1]^Kmax)

#      MH_ratio -= (n_clust + params(obj_prior.prob)[1] - 2.0) * log(obj_mixture_mcmc.prob[1]) + params(obj_prior.prob)[2] * log(1.0 - obj_mixture_mcmc.prob[1]) - log(1.0 - obj_mixture_mcmc.prob[1]^Kmax)

#      if rand(Uniform(0.0,1.0))< exp(MH_ratio)
#        obj_mixture_mcmc.prob[1] = obj_mixture_prop.prob[1]
#      else
        
#      end

#    end
#    obj_mixture_prop.prob[1] = obj_mixture_mcmc.prob[1]
#  end
#  #println("prob", obj_mixture_mcmc.prob[1])
#  for isim in 1:20


#    MH_ratio = 0.0

#    if n_clust == 1
#      n_clust_prop = sample(1:2,1)[1]
#      log_prob_prop = -log(2.0) - log(Kmax - n_clust)
#    elseif n_clust == Kmax
#      n_clust_prop = sample([Kmax, Kmax-1], 1)[1]
#      log_prob_prop = -log(2.0) - log(Kmax - 1.0)
#    else
#      n_clust_prop = sample([n_clust - 1, n_clust, n_clust + 1], 1)[1]
#      log_prob_prop = -log(3.0)
#      if n_clust_prop < n_clust
#        log_prob_prop +=  - log(1.0*n_clust_prop)
#      elseif n_clust_prop > n_clust
#        log_prob_prop += -log(1.0 * (Kmax-n_clust))
#      end
#    end

#    if n_clust_prop == 1
#      log_prob_mcmc = -log(2.0) - log(Kmax - n_clust_prop)
#    elseif n_clust_prop == Kmax
#      log_prob_mcmc = -log(2.0) - log(Kmax - 1.0)
#    else
#      log_prob_mcmc = -log(3.0)
#      if n_clust < n_clust_prop
#        log_prob_mcmc += -log(1.0 * n_clust)
#      elseif n_clust > n_clust_prop
#        log_prob_mcmc += -log(1.0 * (Kmax - n_clust_prop))
#      end
#    end
#    obj_mixture_prop.K[1] = n_clust_prop

#    if n_clust_prop != n_clust
#      MH_ratio = 0.0
#      obj_mixture_prop.K[1] = n_clust_prop
      
#      if n_clust_prop > n_clust

#        is_new = false

#        while is_new == false
#          #println("is_new", is_new)
#          miss_edge_prop[n_clust_prop-1, :] = obj_graph_prop.graph_st[sample(1:size(obj_graph_prop.graph_st, 1))[1], :]
#          #println(miss_edge_prop[n_clust_prop-1, :])
#          #println(miss_edge_prop[1:(n_clust_prop-2), :])
#          is_new = true
#          if n_clust_prop >2
#            for ik in 1:(n_clust_prop-2)

#              if (miss_edge_prop[n_clust_prop-1, 1] == miss_edge_prop[ik, 1]) &  (miss_edge_prop[n_clust_prop-1, 2] == miss_edge_prop[ik, 2])
#                is_new = false
#              end
#              if (miss_edge_prop[n_clust_prop-1, 2] == miss_edge_prop[ik, 1]) & (miss_edge_prop[n_clust_prop-1, 1] == miss_edge_prop[ik, 2])
#                is_new = false
#              end


#            end
#          end
        
#        end

#        it_worked_zeta = update_zeta(obj_mixture_prop, obj_graph_prop)

#        sampling_from_prior(iterations,
#          n_clust_prop,
#          obj_graph_mcmc ,
#          obj_graph_prop ,
#          obj_mixture_mcmc ,
#          obj_mixture_prop ,
#          obj_data_mcmc ,
#          obj_data_prop ,
#          obj_prior)
#      #  println("A")
#        #println([sort!(unique(obj_mixture_mcmc.cluster)), sort!(unique(obj_mixture_prop.cluster))])
#        mapping_mat = matching_cluster_different_sizes!(best_perm, freq_table, obj_mixture_mcmc.cluster, obj_mixture_prop.cluster)
#        map_sep_clust!(mapping_mat, obj_mixture_prop.sep_clust, n_clust)
#        #  println([sort!(unique(obj_mixture_mcmc.cluster)), sort!(unique(obj_mixture_prop.cluster))])
#        for k in 1:n_clust_prop
#          update_which(obj_data_prop, obj_mixture_prop, k)
#          update_param_cluster(obj_data_prop, obj_mixture_prop, k)
#        end

#        for k in 1:n_clust
#          MH_ratio += obj_data_prop.log_likelihood[k] - obj_data_mcmc.log_likelihood[k]
#        end
#        MH_ratio += obj_data_prop.log_likelihood[n_clust_prop]

#      #  println([MH_ratio, (n_clust_prop - 1) * log(obj_mixture_mcmc.prob[1]), (n_clust - 1) * log(obj_mixture_mcmc.prob[1])])
#        # priors
#        MH_ratio += (n_clust_prop - 1) * log(obj_mixture_mcmc.prob[1]) - (n_clust - 1) * log(obj_mixture_mcmc.prob[1])

#        # proposal
#        MH_ratio += log_prob_mcmc - log_prob_prop
      
#        if rand(Uniform(0.0,1.0)) < exp(MH_ratio)
#          do_accept = true
#        else
#          do_accept = false
#        end


#      end

#      if n_clust_prop < n_clust
#        MH_ratio = 0.0
#        if n_clust_prop > 1
#          select_remove = sample(1:(n_clust_prop-1), 1)[1]
#          #println("select_remove", select_remove)
#          if n_clust_prop > 1
#            for imiss in select_remove:(n_clust_prop-1)
#              miss_edge_prop[imiss, :] = miss_edge_mcmc[imiss+1, :]
#            end
#          end
#        end
        
#        it_worked_zeta = update_zeta(obj_mixture_prop, obj_graph_prop)
#        #println("A")
#        #println([sort!(unique(obj_mixture_mcmc.cluster)), sort!(unique(obj_mixture_prop.cluster))])
#        mapping_mat = matching_cluster_different_sizes!(best_perm, freq_table, obj_mixture_mcmc.cluster, obj_mixture_prop.cluster)
#        map_sep_clust!(mapping_mat, obj_mixture_prop.sep_clust, n_clust_prop)


#        mapping_mat = change_label_cluster_and_parameters!(obj_data_prop, obj_mixture_mcmc, obj_mixture_prop.cluster, obj_data_mcmc)
#        map_sep_clust!(mapping_mat, obj_mixture_prop.sep_clust, n_clust_prop)
        

      
#    #    println([sort!(unique(obj_mixture_mcmc.cluster)), sort!(unique(obj_mixture_prop.cluster))])
        
#        for k in 1:n_clust_prop
#          update_which(obj_data_prop, obj_mixture_prop, k)
#          update_param_cluster(obj_data_prop, obj_mixture_prop, k)
#        end

#        for k in 1:(n_clust-1)
#          MH_ratio += obj_data_prop.log_likelihood[k] - obj_data_mcmc.log_likelihood[k]
#        end
#        MH_ratio += -obj_data_mcmc.log_likelihood[n_clust]

#        #println([MH_ratio, (n_clust_prop - 1) * log(obj_mixture_mcmc.prob[1]), (n_clust - 1) * log(obj_mixture_mcmc.prob[1]), log_prob_mcmc, log_prob_prop])
#        # priors

        
#        MH_ratio += (n_clust_prop - 1) * log(obj_mixture_mcmc.prob[1]) - (n_clust - 1) * log(obj_mixture_mcmc.prob[1])

#        # proposal
#        MH_ratio += log_prob_mcmc - log_prob_prop

#        if rand(Uniform(0.0, 1.0)) < exp(MH_ratio)
#          do_accept = true
#        else
#          do_accept = false
#        end


#      end

#    else

#      do_accept = false

#    end

#    if do_accept
#      #  println("accept")
#      for k = 1:obj_mixture_prop.Kmax
#        copy_from_to(obj_data_prop, obj_data_mcmc, k)
#      end
#      copy_from_to(obj_mixture_prop, obj_mixture_mcmc)
#      copy_from_to(obj_graph_prop, obj_graph_mcmc)
      
#    else
  
#      #println("reject")
#      for k = 1:obj_mixture_prop.Kmax
#        copy_from_to(obj_data_mcmc, obj_data_prop, k)
#      end
#      copy_from_to(obj_mixture_mcmc, obj_mixture_prop)
#      copy_from_to(obj_graph_mcmc, obj_graph_prop)

#    end
#    #println("From ", n_clust, " ", n_clust_prop)
#    #println(do_accept)
#    #println((obj_mixture_mcmc.K[1]))
#    #println((miss_edge_mcmc))
#    if (obj_mixture_mcmc.K[1])  <= size(miss_edge_mcmc, 1)
#      for imiss in (obj_mixture_mcmc.K[1]):size(miss_edge_mcmc, 1)
#        miss_edge_mcmc[imiss, :] .= 0
#      end
#    end
#    n_clust = obj_mixture_mcmc.K[1]
#    #println("A")
#    #println(obj_mixture_mcmc.K[1], " ", obj_mixture_prop.K[1])
#    #println(obj_data_mcmc.n_which_obs[1:obj_mixture_mcmc.K[1]], " ", obj_data_prop.n_which_obs[1:obj_mixture_prop.K[1]])
#    #println("B")
#    #println("Post")
#    #println((obj_mixture_mcmc.K[1]))
#    #println((miss_edge_mcmc))
#  end

#end


