abstract type GeneralMixture end

struct TestMixture_V5 <: GeneralMixture

    Kmax::Int64
    n_points::Int64

    cluster::Vector{Int64}
    K::Vector{Int64}
    miss_edge::Matrix{Int64}
    sep_clust::Matrix{Int64}

    root_leaf::Matrix{Int64}
    prob::Vector{Float64}

    cluster_unordered::Vector{Int64}
    vector_unordered::Vector{Int64}

    function TestMixture_V5(
        Kmax::Int64,
        n_points::Int64,
        cluster::Vector{Int64},
        K::Vector{Int64},
        miss_edge::Matrix{Int64},
        sep_clust::Matrix{Int64},
        root_leaf::Matrix{Int64},
        prob::Vector{Float64},
        cluster_unordered::Vector{Int64},
        vector_unordered::Vector{Int64}
    )

        new(Kmax, n_points, cluster, K, miss_edge, sep_clust, root_leaf, prob, cluster_unordered, vector_unordered)

    end

end

function copy_from_to(from::TestMixture_V5, to::TestMixture_V5)

    to.K .= from.K
    to.cluster .= from.cluster
    to.miss_edge .= from.miss_edge
    to.sep_clust .= from.sep_clust
    to.root_leaf .= from.root_leaf
    to.cluster_unordered .= from.cluster_unordered
    to.vector_unordered .= from.vector_unordered
    to.prob .= from.prob
end

function TestMixture_V5(;
    Kmax::Int64,
    miss_edge_init::Matrix{Int64},
    obj_graph::GraphCluter_Vers6,
    K_init::Int64 
)
    n_points = obj_graph.n_points
    K = ones(Int64, 1) * K_init
    miss_edge = zeros(Int64, Kmax-1, 2)
    root_leaf = zeros(Int64, n_points, 2)
    miss_edge[1:size(miss_edge_init, 1), :] = miss_edge_init

    sep_clust = zeros(Int64, Kmax - 1, 2)

    cluster = zeros(Int64, n_points)
    cluster_unordered = zeros(Int64, n_points)
    vector_unordered = zeros(Int64, Kmax)

    out = TestMixture_V5(Kmax, n_points, cluster, K, miss_edge, sep_clust, root_leaf, [0.5], cluster_unordered, vector_unordered)

    it_worked = update_zeta(out, obj_graph)
    if it_worked == false
        println("Wrong Number of clusters")
    end

    return out, it_worked

end

function update_zeta(obj_clust::TestMixture_V5, obj_graph::GraphCluter_Vers6)

    start_index::Int64 = obj_graph.graph_st[1, 1]
    index_k = [1]
    index_iter = [1]
    
    update_zeta_iter(start_index, obj_clust, obj_graph, index_k, index_iter, 1)

    return index_k == obj_clust.K

end



function update_zeta_iter(
    obs::Int64,
    obj_clust::TestMixture_V5,
    obj_graph::GraphCluter_Vers6,
    index_k::Vector{Int64},
    index_iter::Vector{Int64},
    current_k::Int64,
)

    index_iter[1] = index_iter[1] + 1
    if index_iter[1] > 2.0 * (obj_clust.n_points)
        return nothing
    end

    miss_edge = obj_clust.miss_edge
    sep_clust = obj_clust.sep_clust
    zeta = obj_clust.cluster
    neigh_st = obj_graph.neigh_st


    next_k::Int64 = current_k
    if obs == 0
        println("obs", miss_edge)
        println("a = ", obj_graph.n_neigh_st)
        println("b = ", obj_graph.neigh_st)
    end
    zeta[obs] = current_k
    if obj_graph.n_neigh_st[obs] > 0
        for ibranch = 1:obj_graph.n_neigh_st[obs]


            new_node::Int64 = neigh_st[obs][ibranch]

            if obj_clust.K[1] >1
                for imiss = 1:(obj_clust.K[1]-1)

                    if ((miss_edge[imiss, 1] == obs) & (miss_edge[imiss, 2] == new_node)) |
                       ((miss_edge[imiss, 1] == new_node) & (miss_edge[imiss, 2] == obs))

                        sep_clust[imiss, 1] = current_k
                        sep_clust[imiss, 2] = index_k[1] + 1
                        index_k[1] = index_k[1] + 1
                        next_k = index_k[1]
                        #if obj_graph.n_neigh_st[new_node] > 0


                        #end

                    end

                end
            end
            #for imiss = 1:obj_clust.K[1]

            #    if ((miss_edge[imiss, 1] == obs) & (miss_edge[imiss, 2] == new_node)) |
            #       ((miss_edge[imiss, 1] == new_node) & (miss_edge[imiss, 2] == obs))

            #        sep_clust[imiss, 1] = current_k
            #        sep_clust[imiss, 2] = index_k[1] + 1
            #        index_k[1] = index_k[1] + 1
            #        next_k = index_k[1]
            #        #if obj_graph.n_neigh_st[new_node] > 0


            #        #end

            #    end

            #end
            

            update_zeta_iter(new_node, obj_clust, obj_graph, index_k, index_iter, next_k)
            if index_iter[1] > 2.0 * (obj_clust.n_points)
                return nothing
            end
            next_k = current_k
        end
        return nothing
    end


end




function random_matching_cluster!(best_perm::Vector{Int64}, freq_table::Matrix{Int64}, cluster_orig::Vector{Int64}, cluster_to_transform::Vector{Int64})::Dict{Int64,Int64}

    un_orig = sort!(unique(cluster_orig))
    un_prop = sort!(unique(cluster_to_transform))



    #@assert length(un_orig) == length(un_prop)
    #freq_table .= 0
    #for irow = 1:size(un_orig, 1)
    #    for icol = 1:size(un_prop, 1)
    #        freq_table[irow, icol] = sum((cluster_orig .== un_orig[irow]) .& (cluster_to_transform .== un_prop[icol]))
    #    end
    #end

    #best_sum::Float64 = -Inf
    #best_perm .= 0
    #for p in permutations(1:size(un_prop, 1))
    #    s = sum(freq_table[i, p[i]] for i in 1:size(un_orig, 1))
    #    if s > best_sum
    #        best_sum = s
    #        best_perm .= p
    #    end
    #end
    mapping = Dict(zip(un_prop[sample(1:length(un_prop), length(un_prop) , replace=false)], un_orig))
    #mapping = Dict(zip(vec_from, vec_to))
    #println("mapping=", mapping)
    cluster_to_transform .= getindex.(Ref(mapping), cluster_to_transform)

    return mapping

end

function as_perfect_as_possible_matching_cluster!(best_perm::Vector{Int64}, freq_table::Matrix{Int64}, cluster_orig::Vector{Int64}, cluster_to_transform::Vector{Int64})::Dict{Int64,Int64}

    un_orig = sort!(unique(cluster_orig))
    un_prop = sort!(unique(cluster_to_transform))
    @assert length(un_orig) == length(un_prop)
    freq_table .= 0
    for irow = 1:size(un_orig, 1)
        for icol = 1:size(un_prop, 1)
            freq_table[irow, icol] = sum((cluster_orig .== un_orig[irow]) .& (cluster_to_transform .== un_prop[icol]))
        end
    end

    best_sum::Float64 = -Inf
    best_perm .= 0
    for p in permutations(1:size(un_prop, 1))
        s = sum(freq_table[i, p[i]] for i in 1:size(un_orig, 1))
        if s > best_sum
            best_sum = s
            best_perm .= p
        end
    end
    mapping = Dict(zip(un_prop[best_perm], un_orig))
    #mapping = Dict(zip(vec_from, vec_to))
    #println("mapping=", mapping)
    cluster_to_transform .= getindex.(Ref(mapping), cluster_to_transform)

    return mapping

end

function matching_cluster!(best_perm::Vector{Int64}, freq_table::Matrix{Int64}, cluster_orig::Vector{Int64}, cluster_to_transform::Vector{Int64})::Dict{Int64,Int64}
    
    #return as_perfect_as_possible_matching_cluster!(best_perm, freq_table, cluster_orig, cluster_to_transform)

    return random_matching_cluster!(best_perm, freq_table, cluster_orig, cluster_to_transform)
end
    
    

function map_sep_clust!(mapping::Dict{Int64,Int64}, sep_clust::Matrix{Int64}, nclust::Int64)

    for irow = 1:(nclust-1)
        sep_clust[irow, 1] = mapping[sep_clust[irow, 1]]
        sep_clust[irow, 2] = mapping[sep_clust[irow, 2]]
    end

    return nothing

end



function matching_cluster_different_sizes_deterministic!(best_perm::Vector{Int64}, freq_table::Matrix{Int64}, cluster_orig::Vector{Int64}, cluster_to_transform::Vector{Int64})::Dict{Int64,Int64}

    un_orig = sort!(unique(cluster_orig))
    un_prop = sort!(unique(cluster_to_transform))
    
    #if size(un_orig, 1) > size(un_prop,1)

    #    push!(un_prop, 100)

    #end
    #if size(un_orig, 1) < size(un_prop, 1)
    #    push!(un_orig, 100)
    #end
    if size(un_orig, 1) > size(un_prop,1)

        push!(un_prop, size(un_orig, 1))

    end
    if size(un_orig, 1) < size(un_prop, 1)
        push!(un_orig, size(un_prop, 1))
    end
    #println("un_orig=", un_orig)
    #println("un_prop=", un_prop)

    freq_table .= 0
    for irow = 1:size(un_orig, 1)
        for icol = 1:size(un_prop, 1)
            freq_table[irow, icol] = sum((cluster_orig .== un_orig[irow]) .& (cluster_to_transform .== un_prop[icol]))
        end
    end

    best_sum::Float64 = -Inf
    best_perm .= 0
    for p in permutations(1:size(un_prop, 1))
        s = sum(freq_table[i, p[i]] for i in 1:min(size(un_orig, 1), size(un_prop, 1)))
        if s > best_sum
            best_sum = s
            best_perm[1:size(un_prop, 1)] = p
        end
    end
    mapping = Dict(zip(un_prop[best_perm[1:size(un_prop, 1)]], un_orig))
    #mapping = Dict(zip(vec_from, vec_to))
    #println("mapping=", mapping)
    cluster_to_transform .= getindex.(Ref(mapping), cluster_to_transform)

    return mapping

end



function matching_cluster_different_sizes_random!(best_perm::Vector{Int64}, freq_table::Matrix{Int64}, cluster_orig::Vector{Int64}, cluster_to_transform::Vector{Int64})::Dict{Int64,Int64}

    un_orig = sort!(unique(cluster_orig))
    un_prop = sort!(unique(cluster_to_transform))

    #if size(un_orig, 1) > size(un_prop,1)

    #    push!(un_prop, 100)

    #end
    #if size(un_orig, 1) < size(un_prop, 1)
    #    push!(un_orig, 100)
    #end
    if size(un_orig, 1) > size(un_prop, 1)

        push!(un_prop, size(un_orig, 1))

    end
    if size(un_orig, 1) < size(un_prop, 1)
        push!(un_orig, size(un_prop, 1))
    end
    #println("un_orig=", un_orig)
    #println("un_prop=", un_prop)

    #freq_table .= 0
    #for irow = 1:size(un_orig, 1)
    #    for icol = 1:size(un_prop, 1)
    #        freq_table[irow, icol] = sum((cluster_orig .== un_orig[irow]) .& (cluster_to_transform .== un_prop[icol]))
    #    end
    #end

    #best_sum::Float64 = -Inf
    #best_perm .= 0
    #for p in permutations(1:size(un_prop, 1))
    #    s = sum(freq_table[i, p[i]] for i in 1:min(size(un_orig, 1), size(un_prop, 1)))
    #    if s > best_sum
    #        best_sum = s
    #        best_perm[1:size(un_prop, 1)] = p
    #    end
    #end

    #un_prop[sample(1:length(un_prop), length(un_prop), replace=false)]
    mapping = Dict(zip(un_prop[sample(1:length(un_prop), length(un_prop), replace=false)], un_orig))
    #mapping = Dict(zip(vec_from, vec_to))
    #println("mapping=", mapping)
    cluster_to_transform .= getindex.(Ref(mapping), cluster_to_transform)

    return mapping

end

function matching_cluster_different_sizes!(best_perm::Vector{Int64}, freq_table::Matrix{Int64}, cluster_orig::Vector{Int64}, cluster_to_transform::Vector{Int64})::Dict{Int64,Int64}

    return matching_cluster_different_sizes_random!(best_perm, freq_table, cluster_orig, cluster_to_transform)
end
#function update_zeta_from_known_clusters_sep(obj_clust::TestMixture_V5, obj_graph::GraphCluter_Vers6)

#    start_index::Int64 = obj_graph.graph_st[1, 1]
#    index_k = [1]
#    index_iter = [1]

#    update_zeta_iter(start_index, obj_clust, obj_graph, index_k, index_iter, obj_clust.Kmax +1)

#    return index_k == obj_clust.K

#end



#function update_zeta_from_known_clusters_sep_iter(
#    obs::Int64,
#    obj_clust::TestMixture_V5,
#    obj_graph::GraphCluter_Vers6,
#    index_k::Vector{Int64},
#    index_iter::Vector{Int64},
#    current_k::Int64,
#)

#    index_iter[1] = index_iter[1] + 1
#    if index_iter[1] > 2.0 * (obj_clust.n_points)
#        return nothing
#    end

#    miss_edge = obj_clust.miss_edge
#    sep_clust = obj_clust.sep_clust
#    zeta = obj_clust.cluster
#    neigh_st = obj_graph.neigh_st


#    next_k::Int64 = current_k
#    zeta[obs] = current_k
#    if obj_graph.n_neigh_st[obs] > 0
#        for ibranch = 1:obj_graph.n_neigh_st[obs]


#            new_node::Int64 = neigh_st[obs][ibranch]

#            for imiss = 1:obj_clust.K[1]

#                if ((miss_edge[imiss, 1] == obs) & (miss_edge[imiss, 2] == new_node)) |
#                   ((miss_edge[imiss, 1] == new_node) & (miss_edge[imiss, 2] == obs))

#                    sep_clust[imiss, 1] = current_k
#                    sep_clust[imiss, 2] = index_k[1] + 1
#                    index_k[1] = index_k[1] + 1
#                    next_k = index_k[1]
#                    #if obj_graph.n_neigh_st[new_node] > 0


#                    #end

#                end

#            end

#            update_zeta_iter(new_node, obj_clust, obj_graph, index_k, index_iter, next_k)
#            if index_iter[1] > 2.0 * (obj_clust.n_points)
#                return nothing
#            end
#            next_k = current_k
#        end
#        return nothing
#    end


#end


#obj_clust.cluster

#  for iset in 1:obj_graph

#    neigh_st
#  end
