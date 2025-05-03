
function sampling_w(
  iterations::Int64,
  obj_graph_mcmc::GraphCluter_Vers5,
  obj_graph_prop::GraphCluter_Vers5,
  obj_mixture_mcmc::TestMixture_V5,
  obj_mixture_prop::TestMixture_V5,
  obj_data_mcmc::TD,
  obj_data_prop::TD
) where {TD<:GeneralData}

  n_clust::Int64 = obj_mixture_mcmc.K[1]
  icol::Int64 = 0
  sd_vector::Vector{Float64} = [0.0001, 0.001, 0.01, 0.5, 1.0]
  z_sample::Float64 = 0.0
  index_visited::Vector{Int64} = zeros(Int64, obj_graph_mcmc.n_points)
  bool_is_the_same::Bool = false
  iter_app::Int64 = 0
  it_worked_zeta::Bool = false
  do_accept::Bool = false

  node_a::Int64 = 0
  node_b::Int64 = 0

  MH_ratio::Float64 = 0.0

  freq_table::Matrix{Int64} = zeros(Float64, n_clust, n_clust)
  best_sum::Float64 = 0.0
  best_perm::Vector{Int64} = zeros(Float64, n_clust)
  for irow = 1:size(obj_graph_mcmc.neigh_graph, 1)

    for icol_app = 1:size(obj_graph_mcmc.neigh_graph[irow], 1)

      icol = obj_graph_mcmc.neigh_graph[irow][icol_app]

      if irow < icol
        z_sample = sample(sd_vector, 1)[1]
        obj_graph_prop.weight_mat.data[irow, icol] = rand(Normal(obj_graph_mcmc.weight_mat.data[irow, icol], z_sample))
        obj_graph_prop.weight_mat.data[icol, irow] = obj_graph_prop.weight_mat.data[irow, icol]
        if (obj_graph_prop.weight_mat.data[irow, icol] > 0.0) & (obj_graph_prop.weight_mat.data[irow, icol] < 1.0)

          index_visited .= 0
          update_st(obj_graph_prop, index_visited)

          if n_clust > 1
            ## controllo se è compatibile con la partizione esistente
            bool_is_the_same = true

            it_worked_zeta = update_zeta(obj_mixture_prop, obj_graph_prop)
          
            # controllo che i divisori sono possibili
            for ik = 1:(n_clust-1)

              node_a = obj_mixture_prop.miss_edge[ik, 1]
              node_b = obj_mixture_prop.miss_edge[ik, 2]

              bool_is_the_same = bool_is_the_same & ((node_b in obj_graph_prop.neigh_st[node_a]) || (node_a in obj_graph_prop.neigh_st[node_b]))

              #println(bool_is_the_same)

            end
            

    

            if bool_is_the_same
              if it_worked_zeta == false
                println("zz=", obj_mixture_prop.cluster)
                println("a=", obj_mixture_prop.miss_edge)
                println("b=", obj_mixture_prop.sep_clust)
                println("zz=", obj_mixture_mcmc.cluster)
                println("a=", obj_mixture_mcmc.miss_edge)
                println("b=", obj_mixture_mcmc.sep_clust)
                error("didn't work")
              end
              # controllo zeta
              bool_is_the_same = randindex(obj_mixture_prop.cluster, obj_mixture_mcmc.cluster)[2] == 1

              if bool_is_the_same

                obj_mixture_prop.cluster .= obj_mixture_mcmc.cluster
                do_accept = true
              else
                # devo calcolare il rapporto metropolis, ma prima devo uniformare le zeta in qualche modo

                mapping_mat = matching_cluster!(best_perm, freq_table, obj_mixture_mcmc.cluster, obj_mixture_prop.cluster)
                map_sep_clust!(mapping_mat, obj_mixture_prop.sep_clust, n_clust)
                
          #error("")
                
                ## Metropolis

                MH_ratio = 0.0

                for k in 1:n_clust
                  update_which(obj_data_prop, obj_mixture_prop, k)
                  update_param_cluster(obj_data_prop, obj_mixture_prop, k)

                  update_which(obj_data_mcmc, obj_mixture_mcmc, k)
                  update_param_cluster(obj_data_mcmc, obj_mixture_mcmc, k)

                  MH_ratio += obj_data_prop.log_likelihood[k] - obj_data_mcmc.log_likelihood[k]

                end
                #println(MH_ratio)
                
                if rand(Uniform(0.0,1.0)) < exp(MH_ratio)
                  do_accept = true
                else
                  do_accept = false
                end

              end

            else

              do_accept = false

            end



          else

            do_accept = true

          end
        else

          do_accept = false
        end

        ### QUI CAMBIO I PARAMRETRI SE ACCETTO O RIFIUTO
#@time begin
  #println("do_accept=", do_accept)
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
#end
      end

    end


  end




end
