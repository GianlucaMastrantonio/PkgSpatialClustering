function spatial_cluster_v1(;
    iter_mcmc::@NamedTuple{iterations::Int64, burnin::Int64, thin::Int64, n_chains::Int64}=(
        iterations=1000,
        burnin=100,
        thin=1,
        n_chains=1,
    ),
    obj_graph::GraphCluter_Vers5,
    obj_cluster::TestMixture_V5,
    obj_data::TD,
    obj_prior::PriorsMod1_V4,
) where {TD<:GeneralData}

    # indices
    n_points = obj_data.n_points


    obj_graph_mcmc::GraphCluter_Vers5 = GraphCluter_Vers5(graph=obj_graph.graph, weight_mat=deepcopy(obj_graph.weight_mat), coords=obj_graph.coords)
    obj_graph_prop::GraphCluter_Vers5 = GraphCluter_Vers5(graph=obj_graph.graph, weight_mat=deepcopy(obj_graph.weight_mat), coords=obj_graph.coords)

    obj_mixture_mcmc, _ = TestMixture_V5(Kmax=obj_cluster.Kmax, miss_edge_init=deepcopy(obj_cluster.miss_edge), obj_graph=obj_graph_mcmc, K_init = obj_cluster.K[1])

    obj_mixture_prop, _ = TestMixture_V5(Kmax=obj_cluster.Kmax, miss_edge_init=deepcopy(obj_cluster.miss_edge), obj_graph=obj_graph_prop, K_init=obj_cluster.K[1])

    if typeof(obj_data) == GpDataMarginalized_Vers9
        
        obj_data_mcmc = GpDataMarginalized_Vers9(obs=obj_data.obs, gp_init=deepcopy(obj_data.gp), mu_init=deepcopy(obj_data.mu), sigma2_init=deepcopy(obj_data.sigma2), tau2_init=deepcopy(obj_data.tau2), rho_init=deepcopy(obj_data.rho), obj_mixture=obj_mixture_mcmc, obj_graph=obj_graph_mcmc)

        obj_data_prop = GpDataMarginalized_Vers9(obs=obj_data.obs, gp_init=deepcopy(obj_data.gp), mu_init=deepcopy(obj_data.mu), sigma2_init=deepcopy(obj_data.sigma2), tau2_init=deepcopy(obj_data.tau2), rho_init=deepcopy(obj_data.rho), obj_mixture=obj_mixture_mcmc, obj_graph=obj_graph_mcmc)
    else
        obj_data_mcmc = GpData_Vers9(obs=obj_data.obs, gp_init=deepcopy(obj_data.gp), mu_init=deepcopy(obj_data.mu), sigma2_init=deepcopy(obj_data.sigma2), tau2_init=deepcopy(obj_data.tau2), rho_init=deepcopy(obj_data.rho), obj_mixture=obj_mixture_mcmc, obj_graph=obj_graph_mcmc)

        obj_data_prop = GpData_Vers9(obs=obj_data.obs, gp_init=deepcopy(obj_data.gp), mu_init=deepcopy(obj_data.mu), sigma2_init=deepcopy(obj_data.sigma2), tau2_init=deepcopy(obj_data.tau2), rho_init=deepcopy(obj_data.rho), obj_mixture=obj_mixture_mcmc, obj_graph=obj_graph_mcmc)
    end
    

    #obj_data_mcmc::GpData_Vers9 = GpData_Vers9(obs=obj_data.obs, gp_init=deepcopy(obj_data.gp), mu_init=deepcopy(obj_data.mu), sigma2_init=deepcopy(obj_data.sigma2), tau2_init=deepcopy(obj_data.tau2), rho_init=deepcopy(obj_data.rho), obj_mixture=obj_mixture_mcmc, obj_graph=obj_graph_mcmc)

    #obj_data_prop = GpData_Vers9(obs=obj_data.obs, gp_init=deepcopy(obj_data.gp), mu_init=deepcopy(obj_data.mu), sigma2_init=deepcopy(obj_data.sigma2), tau2_init=deepcopy(obj_data.tau2), rho_init=deepcopy(obj_data.rho), obj_mixture=obj_mixture_prop, obj_graph=obj_graph_prop)

    
    ## iterations
    iterations = iter_mcmc.iterations
    burnin = iter_mcmc.burnin
    thin = iter_mcmc.thin
    sample_to_save = Int64(floor((iterations - burnin) / thin))
    thin_burnin = burnin
    iterations::Int64 = 0
    p2 = Progress(
        burnin + (sample_to_save - 1) * thin,
        desc="iterations ",
        offset=0,
        showspeed=true,
    )
    
    #p3 = Progress(
    #    burnin + (sample_to_save - 1) * thin,
    #    desc="clusters ",
    #    offset=0,
    #    showspeed=true,
    #)

    # output mcmc
    zeta_out = zeros(Int64, n_points, sample_to_save)
    K_out = zeros(Int64, sample_to_save)
    miss_edge_out = zeros(Int64, sample_to_save, 2 * (obj_mixture_mcmc.Kmax - 1))
    spanning_tree_out = zeros(Int64, n_points - 1, 2, sample_to_save)

    index_node_weight_mat::Vector{Int64} = findall(x -> x != 0, obj_graph_mcmc.weight_mat[:])
    w_index_out = zeros(Float64, size(index_node_weight_mat,1), sample_to_save)

    sep_clust_out = zeros(Int64, obj_mixture_mcmc.Kmax - 1, 2, sample_to_save)
    miss_edge_out = zeros(Int64, obj_mixture_mcmc.Kmax - 1, 2, sample_to_save)
        
    n_cluster_out = zeros(Int64, sample_to_save)

    mu_out = zeros(Float64, obj_mixture_mcmc.Kmax, sample_to_save)
    tau2_out = zeros(Float64, obj_mixture_mcmc.Kmax, sample_to_save)
    sigma2_out = zeros(Float64, obj_mixture_mcmc.Kmax, sample_to_save)
    rho_out = zeros(Float64, obj_mixture_mcmc.Kmax, sample_to_save)
    gp_out = zeros(Float64, n_points, obj_mixture_mcmc.Kmax, sample_to_save)

    for k in 1:obj_mixture_mcmc.K[1]
        #println("k", k)
        update_param_cluster(obj_data_mcmc, obj_mixture_mcmc, k::Int64)
    end

    ### adapt 
    


    for iMCMC = 1:sample_to_save
        for jMCMC = 1:thin_burnin

            iterations += 1
            #println("K = ", obj_mixture_mcmc.K[1], " - ", countmap(obj_mixture_mcmc.cluster))
            #println(countmap(obj_mixture_mcmc.cluster))
            #ProgressMeter.next!(p2; showvalues=[(:iterations, iterations)])
            ProgressMeter.next!(p2; showvalues=[
                (:iterations, iterations),
                (:K, obj_mixture_mcmc.K[1]),
                (:clusters, countmap(obj_mixture_mcmc.cluster))
            ])

            
            #println("B")
            sampling_separator(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop)
            
            #println("A")
            #@time s
            sampling_separator_jump(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop)
            
            
            @toggled_assert sum(obj_data_mcmc.n_which_obs[1:obj_mixture_mcmc.K[1]]) == obj_data_mcmc.n_points
            @toggled_assert sum(obj_data_mcmc.n_which_obs[1:obj_mixture_mcmc.K[1]] .> 0) == obj_mixture_mcmc.K[1]
            @toggled_assert length(sort(vcat([obj_data_mcmc.which_obs[i][1:obj_data_mcmc.n_which_obs[i]] for i in 1:obj_mixture_mcmc.K[1]]...))) == obj_data_mcmc.n_points


            @toggled_assert sum(obj_data_prop.n_which_obs[1:obj_mixture_prop.K[1]]) == obj_data_prop.n_points
            @toggled_assert sum(obj_data_prop.n_which_obs[1:obj_mixture_prop.K[1]] .> 0) == obj_mixture_prop.K[1]
            @toggled_assert length(sort(vcat([obj_data_prop.which_obs[i][1:obj_data_prop.n_which_obs[i]] for i in 1:obj_mixture_prop.K[1]]...))) == obj_data_prop.n_points

            

            
            #println("B")
            #@time 
            sampling_ncluster(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)

            @toggled_assert sum(obj_data_mcmc.n_which_obs[1:obj_mixture_mcmc.K[1]]) == obj_data_mcmc.n_points
            @toggled_assert sum(obj_data_mcmc.n_which_obs[1:obj_mixture_mcmc.K[1]] .> 0) == obj_mixture_mcmc.K[1]
            @toggled_assert length(sort(vcat([obj_data_mcmc.which_obs[i][1:obj_data_mcmc.n_which_obs[i]] for i in 1:obj_mixture_mcmc.K[1]]...))) == obj_data_mcmc.n_points

            #println(obj_mixture_mcmc.K[1], " ",obj_mixture_prop.K[1])
            #println(obj_data_mcmc.n_which_obs[1:obj_mixture_mcmc.K[1]], " ", obj_data_prop.n_which_obs[1:obj_mixture_prop.K[1]])
            @toggled_assert sum(obj_data_prop.n_which_obs[1:obj_mixture_prop.K[1]]) == obj_data_prop.n_points
            @toggled_assert sum(obj_data_prop.n_which_obs[1:obj_mixture_prop.K[1]] .> 0) == obj_mixture_prop.K[1]
            @toggled_assert length(sort(vcat([obj_data_prop.which_obs[i][1:obj_data_prop.n_which_obs[i]] for i in 1:obj_mixture_prop.K[1]]...))) == obj_data_prop.n_points


            

            # ! LASCIARLO QUI
            from_marg_to_full(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
            
            #println("C")
            #@time 
            sampling_w(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
            
            sampling_mu(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
            @toggled_assert obj_data_prop.log_det[1] == obj_data_mcmc.log_det[1]
            
            
            
            sampling_tau2(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
            @toggled_assert obj_data_prop.log_det[1] == obj_data_mcmc.log_det[1]
            
            
            
            
            sampling_sigma2(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
            @toggled_assert obj_data_prop.log_det[1] == obj_data_mcmc.log_det[1]


            
            sampling_rho_and_gp(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
            @toggled_assert obj_data_prop.log_det[1] == obj_data_mcmc.log_det[1]


            # ! LASCIARLO QUI
            from_full_to_marg(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior)
            
            #


            for k in 1:obj_mixture_mcmc.K[1]
                #println("k", k)
                update_param_cluster(obj_data_mcmc, obj_mixture_mcmc, k)
            end

        end
        thin_burnin = thin
        ## sample_to_save
        zeta_out[:, iMCMC] = obj_mixture_mcmc.cluster
        spanning_tree_out[:, :, iMCMC] = obj_graph_mcmc.graph_st        
        w_index_out[:, iMCMC] = obj_graph_mcmc.weight_mat[index_node_weight_mat]
        sep_clust_out[:,:,iMCMC] = obj_mixture_mcmc.sep_clust
        miss_edge_out[:,:,iMCMC] = obj_mixture_mcmc.miss_edge
        n_cluster_out[iMCMC] = obj_mixture_mcmc.K[1]

        mu_out[:, iMCMC] = obj_data_mcmc.mu
        tau2_out[:, iMCMC] = obj_data_mcmc.tau2
        sigma2_out[:, iMCMC] = obj_data_mcmc.sigma2
        rho_out[:, iMCMC] = obj_data_mcmc.rho
        gp_out[:, :, iMCMC] = obj_data_mcmc.gp

        #K_out = zeros(Int64, sample_to_save)
        #miss_edge_out = zeros(Int64, sample_to_save, 2 * (Kmax - 1))
    end

    return (
        iter_mcmc=iter_mcmc,
        index_node_weight_mat = index_node_weight_mat,
        zeta_out=zeta_out,
        spanning_tree_out=spanning_tree_out,
        w_index_out=w_index_out,
        sep_clust_out = sep_clust_out,
        miss_edge_out=miss_edge_out,
        n_cluster_out=n_cluster_out,
        mu_out=mu_out,
        tau2_out=tau2_out,
        sigma2_out=sigma2_out,
        rho_out=rho_out,
        gp_out=gp_out
    )


end
