function sampling_st(
    iterations::Int64,
    obj_graph_mcmc::GraphCluter_Vers6,
    obj_mixture_mcmc::TestMixture_V5,
    k_mcmc::Vector{Int64},
)



end
#function sampling_st(iterations::Int64, obj_graph_mcmc::GraphCluter_Vers6, obj_mixture_mcmc::TestMixture_V5, k_mcmc::Vector{Int64})


#  n_points = obj_graph_mcmc.n_points
#  K = k_mcmc[1]
#  index_visited = zeros(Int64, n_points)
#  if K >1

#    mixture_graphs_st = Vector{SimpleGraph{Int64}}(undef, K)
#    mixture_roots = zeros(Int64, K)
#    mixture_index = Vector{Vector{Int64}}(undef, K)

#    app_graph_st = SimpleGraph(n_points)
#    for k in 1:K
#      println("k =", k)
#      mixture_index[k] = findall(x -> x == k, obj_mixture_mcmc.cluster)
#      println(size(mixture_index[k], 1))

#      sub_graph = induced_subgraph(obj_graph_mcmc.graph, mixture_index[k])

#      #if iterations == 74
#      #  println(dump(sub_graph))
#      #  println(is_connected(sub_graph[1]))
#      #  mixture_graphs_st, mixture_roots[k] = wilsons_algorithm_weighted_test(sub_graph[1], obj_graph_mcmc.weight_mat.data[mixture_index[k], mixture_index[k]])
#      #end
#      mixture_graphs_st, mixture_roots[k] = wilsons_algorithm_weighted(sub_graph[1], obj_graph_mcmc.weight_mat.data[mixture_index[k], mixture_index[k]])

#      for e = 1:size(mixture_graphs_st.fadjlist, 1)

#        app_graph_st.fadjlist[mixture_index[k][e]] = mixture_index[k][mixture_graphs_st.fadjlist[e]]

#      end

#    end
#    app_graph_st.ne = n_points - 1
#    for k in 1:(K-1)
#      app_graph_st.fadjlist[obj_mixture_mcmc.miss_edge[k, 1]] = sort([app_graph_st.fadjlist[obj_mixture_mcmc.miss_edge[k, 1]]; obj_mixture_mcmc.miss_edge[k, 2]])

#      app_graph_st.fadjlist[obj_mixture_mcmc.miss_edge[k, 2]] = sort([app_graph_st.fadjlist[obj_mixture_mcmc.miss_edge[k, 2]]; obj_mixture_mcmc.miss_edge[k, 1]])
#    end


#    # ordino i punti
#    order_wilsons(app_graph_st, mixture_roots[1], obj_graph_mcmc.graph_st, index_visited)


#    obj_graph_mcmc.n_neigh_st[:] .= 0

#    for ig = 1:(n_points-1)
#      a = obj_graph_mcmc.graph_st[ig, 1]
#      b = obj_graph_mcmc.graph_st[ig, 2]

#      obj_graph_mcmc.n_neigh_st[a] += 1
#      obj_graph_mcmc.neigh_st[a][obj_graph_mcmc.n_neigh_st[a]] = b



#    end

#    ### test

#    save_zeta = deepcopy(obj_mixture_mcmc.cluster)


#    it_worked = update_zeta(obj_mixture_mcmc, obj_graph_mcmc)
#    if it_worked == false
#      println("Wrong Number of clusters")
#    end

#    for k in 1:K
#      obj_mixture_mcmc.cluster[mixture_index[k]] .= save_zeta[mixture_index[k][1]]
#    end
#    #if obj_mixture_mcmc.cluster != save_zeta

#    #  println("a=",  save_zeta)
#    #  println("b=", obj_mixture_mcmc.cluster)
#    #  error("")
#    #end

#    #  println(obj_mixture_mcmc.cluster .== save_zeta)

#    #end
#    #println(obj_mixture_mcmc.cluster == save_zeta)
#    @assert obj_mixture_mcmc.cluster == save_zeta "zeta"
#  else
#    app_graph_st, root = wilsons_algorithm_weighted(obj_graph_mcmc.graph, obj_graph_mcmc.weight_mat.data)
#    order_wilsons(app_graph_st, root, obj_graph_mcmc.graph_st, index_visited)


#    obj_graph_mcmc.n_neigh_st[:] .= 0

#    for ig = 1:(n_points-1)
#      a = obj_graph_mcmc.graph_st[ig, 1]
#      b = obj_graph_mcmc.graph_st[ig, 2]

#      obj_graph_mcmc.n_neigh_st[a] += 1
#      obj_graph_mcmc.neigh_st[a][obj_graph_mcmc.n_neigh_st[a]] = b



#    end
#  end



#end
