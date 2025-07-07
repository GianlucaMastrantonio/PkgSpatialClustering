abstract type GeneralGraph end

struct GraphCluter_Vers6 <: GeneralGraph

    n_points::Int64

    graph::Graphs.SimpleGraphs.SimpleGraph{Int64}
    graph_st::Matrix{Int64}

    weight_mat::Symmetric{Float64,Matrix{Float64}}
    coords::Matrix{Float64}

    neigh_graph::Vector{Vector{Int64}}
    neigh_st::Vector{Vector{Int64}}

    n_neigh_graph::Vector{Int64}
    n_neigh_st::Vector{Int64}
    covariates::Array{Float64,3}
    predictors::Vector{Float64}




    function GraphCluter_Vers6(
        n_points::Int64,
        graph::Graphs.SimpleGraphs.SimpleGraph{Int64},
        graph_st::Matrix{Int64},
        weight_mat::Symmetric{Float64,Matrix{Float64}},
        coords::Matrix{Float64},
        neigh_graph::Vector{Vector{Int64}},
        neigh_st::Vector{Vector{Int64}},
        n_neigh_graph::Vector{Int64},
        n_neigh_st::Vector{Int64},
        covariates::Array{Float64,3},
        predictors::Vector{Float64}
    )

        new(
            n_points,
            graph,
            graph_st,
            weight_mat,
            coords,
            neigh_graph,
            neigh_st,
            n_neigh_graph,
            n_neigh_st,
            covariates,
            predictors
        )

    end

end

function copy_from_to(from::GraphCluter_Vers6, to::GraphCluter_Vers6)

    to.graph_st .= from.graph_st
    to.weight_mat.data .= from.weight_mat.data
    for i = 1:size(from.neigh_st,1)
        to.neigh_st[i] .= from.neigh_st[i]
    end    
    to.n_neigh_st .= from.n_neigh_st
    to.predictors .= from.predictors

    


end


function update_st(graph_obj::GraphCluter_Vers6, index_visited::Vector{Int64})

    index_visited .= 0.0
    graph_obj.graph_st .= 0
    app_graph_st = prim_mst(graph_obj.graph, graph_obj.weight_mat.data)
    order_mst(app_graph_st, graph_obj.graph_st, index_visited)

    for i = 1:size(graph_obj.neigh_st, 1)
        graph_obj.neigh_st[i] .= 0
    end

    a::Int64 = 0
    b::Int64 = 0
    graph_obj.n_neigh_st .= 0
    for ig = 1:(graph_obj.n_points-1)

        a = graph_obj.graph_st[ig, 1]
        b = graph_obj.graph_st[ig, 2]
        #if a   == 0
        #    error("a=0")
        #end
        #if b == 0
        #    error("b=0")
        #end

        graph_obj.n_neigh_st[a] += 1
        graph_obj.neigh_st[a][graph_obj.n_neigh_st[a]] = b

    end
end

function GraphCluter_Vers6(;
    graph::Graphs.SimpleGraphs.SimpleGraph{Int64},
    weight_mat::Symmetric{Float64,Matrix{Float64}},
    coords::Matrix{Float64},
    covariates::Array{Float64,3},
    predictors::Vector{Float64}
)
    
    for i = 2:size(covariates, 2)
        for j in 1:(i-1)
            for k in 1:size(covariates,1)
                covariates[k,i,j] = covariates[k,j,i]     
            end
        end
    end
    n_points = size(weight_mat, 1)
    if size(covariates, 2) != size(covariates, 3)
        error("")
    end
    if size(covariates, 2) != size(coords, 1)
        error("covariates and coords must have the same number of columns")
    end
    #app_graph_st = Graphs.prim_mst(graph, weight_mat)
    #wilsons_algorithm(graph)
    graph_st = zeros(Int64, n_points - 1, 2)
    # campiono il grafo
    #app_graph_st, root = wilsons_algorithm_weighted(graph, weight_mat.data)
    app_graph_st = prim_mst(graph, weight_mat.data)
    #kruskal_mst
    #prim_mst
    index_visited = zeros(Int64, n_points)
    order_mst(app_graph_st, graph_st, index_visited)
    #println(graph_st)
    #for i in 1:(n_points-1)
    #  graph_st[i, 1] = app_graph_st[i].src
    #  graph_st[i, 2] = app_graph_st[i].dst
    #end
    #index_visited = zeros(Int64, n_points)
    ## ordino i punti
    #

    neigh_graph::Vector{Vector{Int64}} = graph.fadjlist
    n_neigh_graph::Vector{Int64} = [size(neigh_graph[i], 1) for i = 1:n_points]


    neigh_st::Vector{Vector{Int64}} =
        [zeros(Int64, size(neigh_graph[i], 1)) for i = 1:n_points]
    n_neigh_st::Vector{Int64} = zeros(Int64, n_points)

    for ig = 1:(n_points-1)
        a = graph_st[ig, 1]
        b = graph_st[ig, 2]

        n_neigh_st[a] += 1
        neigh_st[a][n_neigh_st[a]] = b



    end

    return GraphCluter_Vers6(
        n_points,
        graph,
        graph_st,
        weight_mat,
        coords,
        neigh_graph,
        neigh_st,
        n_neigh_graph,
        n_neigh_st,
        covariates,
        predictors
    )
end


function order_mst(
    g::Array{Graphs.SimpleGraphs.SimpleEdge{Int64},1},
    graph_st::Matrix{Int64},
    index_visited::Vector{Int64},
)

    #root::Int64 = sample(1:size(graph_st, 1),1)[1]
    root::Int64 = g[size(graph_st, 1)].dst
    g_new = SimpleGraph(size(graph_st, 1) + 1)
    for e in g
        add_edge!(g_new, e.src, e.dst)
    end

    index_visited .= 0
    index_pos = [1]
    index_visited[root] = 1
    order_mst_recursive(g_new, root, graph_st, index_visited, index_pos)

end


function order_mst_recursive(
    g::SimpleGraph{Int64},
    root::Int64,
    graph_st::Matrix{Int64},
    index_visited::Vector{Int64},
    index_pos::Vector{Int64},
)


    for ig = 1:size(g.fadjlist[root], 1)

        if index_visited[g.fadjlist[root][ig]] == 0

            graph_st[index_pos[1], 1] = root
            graph_st[index_pos[1], 2] = g.fadjlist[root][ig]

            index_pos[1] += 1
            index_visited[g.fadjlist[root][ig]] += 1

            order_mst_recursive(g, g.fadjlist[root][ig], graph_st, index_visited, index_pos)

        else



        end

    end

    return nothing


end

#function GraphCluter_Vers6(graph::Graphs.SimpleGraphs.SimpleGraph{Int64}, weight_mat::Symmetric{Float64,Matrix{Float64}}, coords::Matrix{Float64})

#  n_points = size(weight_mat, 1)

#  #app_graph_st = Graphs.prim_mst(graph, weight_mat)
#  #wilsons_algorithm(graph)
#  graph_st = zeros(Int64, n_points - 1, 2)
#  # campiono il grafo
#  #app_graph_st, root = wilsons_algorithm_weighted(graph, weight_mat.data)
#  app_graph_st = prim_mst(graph, weight_mat.data)

#  index_visited = zeros(Int64, n_points)
#  # ordino i punti
#  order_wilsons(app_graph_st, root, graph_st, index_visited)

#  neigh_graph::Vector{Vector{Int64}} = graph.fadjlist
#  n_neigh_graph::Vector{Int64} = [size(neigh_graph[i], 1) for i in 1:n_points]

#  neigh_st::Vector{Vector{Int64}} = [zeros(Int64, size(neigh_graph[i], 1)) for i in 1:n_points]
#  n_neigh_st::Vector{Int64} = zeros(Int64, n_points)

#  for ig = 1:(n_points-1)
#    a = graph_st[ig,1]
#    b = graph_st[ig, 2]

#    n_neigh_st[a] += 1
#    neigh_st[a][n_neigh_st[a]] = b



#  end

#  return GraphCluter_Vers6(n_points, graph, graph_st, weight_mat, coords, neigh_graph, neigh_st, n_neigh_graph, n_neigh_st)
#end


function order_wilsons(
    g::SimpleGraph{Int64},
    root::Int64,
    graph_st::Matrix{Int64},
    index_visited::Vector{Int64},
)

    index_visited .= 0
    index_pos = [1]
    index_visited[root] = 1
    order_wilsons_recursive(g, root, graph_st, index_visited, index_pos)

end


function order_wilsons_recursive(
    g::SimpleGraph{Int64},
    root::Int64,
    graph_st::Matrix{Int64},
    index_visited::Vector{Int64},
    index_pos::Vector{Int64},
)


    for ig = 1:size(g.fadjlist[root], 1)

        if index_visited[g.fadjlist[root][ig]] == 0

            graph_st[index_pos[1], 1] = root
            graph_st[index_pos[1], 2] = g.fadjlist[root][ig]

            index_pos[1] += 1
            index_visited[g.fadjlist[root][ig]] += 1

            order_wilsons_recursive(
                g,
                g.fadjlist[root][ig],
                graph_st,
                index_visited,
                index_pos,
            )

        else



        end

    end

    return nothing


end

function wilsons_algorithm(g::SimpleGraph)
    n = nv(g)
    tree = SimpleGraph(n)             # Start an empty graph with n vertices
    in_tree = falses(n)
    root = rand(1:n)
    in_tree[root] = true

    while !all(in_tree)
        start = rand(findall(!, in_tree))
        path = Dict{Int,Int}()
        current = start

        # Random walk until hitting the tree
        while !in_tree[current]
            nbrs = collect(neighbors(g, current))
            next = rand(nbrs)
            path[current] = next
            current = next
        end

        # Retrace and add path edges to the tree
        current = start
        while !in_tree[current]
            next = path[current]
            add_edge!(tree, current, next)
            in_tree[current] = true
            current = next
        end
    end

    return tree, root
end




#function wilsons_algorithm_weighted(g::SimpleGraph, W::Matrix{Float64})
#  n = nv(g)
#  tree = SimpleGraph(n)
#  in_tree = falses(n)
#  root = rand(1:n)
#  in_tree[root] = true

#  while !all(in_tree)
#    start = rand(findall(!, in_tree))
#    path = Dict{Int,Int}()
#    current = start

#    while !in_tree[current]
#      nbrs = collect(neighbors(g, current))
#      weights = [W[current, nb] for nb in nbrs]
#      next = sample(nbrs, Weights(weights))
#      path[current] = next
#      current = next
#    end

#    current = start
#    while !in_tree[current]
#      next = path[current]
#      add_edge!(tree, current, next)
#      in_tree[current] = true
#      current = next
#    end
#  end

#  return tree, root
#end
function wilsons_algorithm_weighted(g::SimpleGraph, W::Matrix{Float64})
    iii2 = 0
    n = nv(g)
    tree = SimpleGraph(n)
    in_tree = falses(n)
    root = rand(1:n)
    in_tree[root] = true
    #println("root = ", root)
    visited = falses(n)
    while !all(in_tree)

        start = rand(findall(!, in_tree))
        current = start
        path = Int[]
        came_from = fill(-1, n)  # Use an array instead of Dict for performance
        #println("ssstart = ", start)
        while !in_tree[current]
            nbrs = collect(neighbors(g, current))
            weights = [W[current, nb] for nb in nbrs]
            #println(nbrs)
            #println(weights)

            next = sample(nbrs, Weights(weights))
            came_from[current] = next
            current = next
            iii2 += 1
            if iii2 > 100003
                println("root", root)
                nbrs = collect(neighbors(g, root))
                weights = [W[root, nb] for nb in nbrs]
                println(nbrs)
                println(weights)
                println("start", start)
                nbrs = collect(neighbors(g, start))
                weights = [W[start, nb] for nb in nbrs]
                println(nbrs)
                println(weights)

                path = dijkstra_shortest_paths(g, root)
                println("path length= ", path.dists[start])
                println()
                nnn = zeros(Int64, 10)
                aa = start
                for iiiii = 1:(path.dists[start])
                    nnn[iiiii] = path.parents[aa][1]
                    aa = nnn[iiiii]
                end
                print(nnn)


                error("too many iterations")


            end
        end

        # Loop-erased path construction
        visited .= false  # Use a boolean array instead of Dict
        node = start
        while !in_tree[node]
            push!(path, node)
            visited[node] = true
            node = came_from[node]
            if visited[node]
                # Loop detected: erase
                idx = findfirst(x -> x == node, path)
                path = path[1:idx-1]
                visited[path] .= true
                node = path[end]
            end
        end

        # Add the loop-erased path to the tree
        node = start
        while !in_tree[node]
            next = came_from[node]
            add_edge!(tree, node, next)
            in_tree[node] = true
            node = next
        end
    end

    return tree, root
end


#function wilsons_algorithm_weighted_test(g::SimpleGraph, W::Matrix{Float64})
#  iii1 =0
#  iii2 = 0
#  iii3 = 0
#  iii4 = 0

#  n = nv(g)
#  tree = SimpleGraph(n)
#  in_tree = falses(n)
#  root = rand(1:n)
#  in_tree[root] = true
#  println(root)
#  #error("sss")
#  visited = falses(n)  # Use a boolean array instead of Dict
#  while !all(in_tree)

#    iii1 += 1
#    println("A1w ", iii1)
#    start = rand(findall(!, in_tree))
#    println("start ", start)

#    current = start
#    path = Int[]
#    came_from = fill(-1, n)  # Use an array instead of Dict for performance

#    while !in_tree[current]

#      iii2 += 1
#      println("A2 ", iii2)
#      nbrs = collect(neighbors(g, current))
#      weights = [W[current, nb] for nb in nbrs]
#      next = sample(nbrs, Weights(weights))
#      came_from[current] = next
#      current = next
#      println(current)
#      if iii2 >10003
#        println(weights)
#        error("")
#      end
#    end
#    error("err")
#    # Loop-erased path construction
#    visited .= false
#    node = start
#    while !in_tree[node]
#      iii3 += 1
#      println("A3 ", iii3)
#      push!(path, node)
#      visited[node] = true
#      node = came_from[node]
#      if visited[node]
#        # Loop detected: erase
#        idx = findfirst(x -> x == node, path)
#        path = path[1:idx-1]
#        visited[path] .= true
#        node = path[end]
#      end
#    end

#    # Add the loop-erased path to the tree
#    node = start
#    while !in_tree[node]
#      iii4 += 1
#      println("A4 ", iii4)
#      next = came_from[node]
#      add_edge!(tree, node, next)
#      in_tree[node] = true
#      node = next
#    end
#  end
#  println("end wilsons")
#  return tree, root
#end
