function spatial_cluster_v1(;
    iter_mcmc::@NamedTuple{iterations::Int64, burnin::Int64, thin::Int64, n_chains::Int64}=(
        iterations=1000,
        burnin=100,
        thin=1,
        n_chains=1
    ),
    obj_graph::GraphCluter_Vers6,
    obj_cluster::TestMixture_V5,
    obj_data::TD,
    obj_prior::PriorsMod1_V6,
    temperatures::Vector{Float64}=[1.0],
    amcmc::Vector{MultiType1_V3},
    cohesion::TC
) where {TD<:GeneralData,TC<:CohesionFunction}

    # indices
    n_points = obj_data.n_points
    n_temps = size(temperatures, 1)


    obj_graph_mcmc::Vector{GraphCluter_Vers6} = [GraphCluter_Vers6(graph=obj_graph.graph, weight_mat=deepcopy(obj_graph.weight_mat), coords=obj_graph.coords, covariates=obj_graph.covariates, predictors=deepcopy(obj_graph.predictors)) for _ in 1:n_temps]
    obj_graph_prop::Vector{GraphCluter_Vers6} = [GraphCluter_Vers6(graph=obj_graph.graph, weight_mat=deepcopy(obj_graph.weight_mat), coords=obj_graph.coords, covariates=obj_graph.covariates, predictors=deepcopy(obj_graph.predictors)) for _ in 1:n_temps]

    obj_mixture_mcmc::Vector{TestMixture_V5} = [TestMixture_V5(Kmax=obj_cluster.Kmax, miss_edge_init=deepcopy(obj_cluster.miss_edge), obj_graph=obj_graph_mcmc[it], K_init=obj_cluster.K[1])[1] for it in 1:n_temps]

    obj_mixture_prop::Vector{TestMixture_V5} = [TestMixture_V5(Kmax=obj_cluster.Kmax, miss_edge_init=deepcopy(obj_cluster.miss_edge), obj_graph=obj_graph_prop[it], K_init=obj_cluster.K[1])[1] for it in 1:n_temps]

    obj_data_mcmc::Vector{GpDataMarginalized_Vers9} = [GpDataMarginalized_Vers9(obs=obj_data.obs, gp_init=deepcopy(obj_data.gp), mu_init=deepcopy(obj_data.mu), sigma2_init=deepcopy(obj_data.sigma2), tau2_init=deepcopy(obj_data.tau2), rho_init=deepcopy(obj_data.rho), obj_mixture=obj_mixture_mcmc[it], obj_graph=obj_graph_mcmc[it], temperature=temperatures[it]) for it in 1:n_temps]

    obj_data_prop::Vector{GpDataMarginalized_Vers9} = [GpDataMarginalized_Vers9(obs=obj_data.obs, gp_init=deepcopy(obj_data.gp), mu_init=deepcopy(obj_data.mu), sigma2_init=deepcopy(obj_data.sigma2), tau2_init=deepcopy(obj_data.tau2), rho_init=deepcopy(obj_data.rho), obj_mixture=obj_mixture_mcmc[it], obj_graph=obj_graph_mcmc[it], temperature=temperatures[it]) for it in 1:n_temps]

    adaptive_mcmc::Vector{Vector{MultiType1_V3}} = [[MultiType1_V3(parameters=amcmc[k].par_acc, cov_init=amcmc[k].cov_init, lambda_init=amcmc[k].lambda_init, diag_mat=amcmc[k].diag_mat, target_alpha=amcmc[k].target_alpha, iter_start=1, iter_end=deepcopy(amcmc[k].iter_end), par_a=deepcopy(amcmc[k].par_a), par_b=deepcopy(amcmc[k].par_b)) for k in 1:size(amcmc, 1)] for it in 1:n_temps]

    obj_cohesion_mcmc = [create_same_object(cohesion) for it in 1:n_temps]
    obj_cohesion_prop = [create_same_object(cohesion) for it in 1:n_temps]

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



    # output mcmc
    zeta_out = zeros(Int64, n_points, sample_to_save)
    K_out = zeros(Int64, sample_to_save)
    miss_edge_out = zeros(Int64, sample_to_save, 2 * (obj_mixture_mcmc[1].Kmax - 1))
    spanning_tree_out = zeros(Int64, n_points - 1, 2, sample_to_save)

    index_node_weight_mat::Vector{Int64} = findall(x -> x != 0, obj_graph_mcmc[1].weight_mat[:])
    w_index_out = zeros(Float64, size(index_node_weight_mat, 1), sample_to_save)

    sep_clust_out = zeros(Int64, obj_mixture_mcmc[1].Kmax - 1, 2, sample_to_save)
    miss_edge_out = zeros(Int64, obj_mixture_mcmc[1].Kmax - 1, 2, sample_to_save)

    n_cluster_out = zeros(Int64, sample_to_save)

    mu_out = zeros(Float64, obj_mixture_mcmc[1].Kmax, sample_to_save)
    tau2_out = zeros(Float64, obj_mixture_mcmc[1].Kmax, sample_to_save)
    sigma2_out = zeros(Float64, obj_mixture_mcmc[1].Kmax, sample_to_save)
    rho_out = zeros(Float64, obj_mixture_mcmc[1].Kmax, sample_to_save)
    gp_out = zeros(Float64, n_points, obj_mixture_mcmc[1].Kmax, sample_to_save)

    predictors_out = zeros(Float64, size(obj_graph_mcmc[1].predictors, 1), sample_to_save)

    for it in 1:size(temperatures, 1)
        for k in 1:obj_mixture_mcmc[it].K[1]
            #println("k", k)
            update_param_cluster(obj_data_mcmc[it], obj_mixture_mcmc[it], k::Int64, temperatures[it])
        end
    end


    ### adapt 




    which_data_to_temp::Vector{Int64} = collect(1:n_temps) ## il primo elemento mi dice quale tra le temperature è associata al dataset 1
    which_first::Int64 = 1
    for iMCMC = 1:sample_to_save
        for jMCMC = 1:thin_burnin

            iterations += 1
            #println(iterations)
            #println("K = ", obj_mixture_mcmc.K[1], " - ", countmap(obj_mixture_mcmc.cluster))
            #println(countmap(obj_mixture_mcmc.cluster))
            #ProgressMeter.next!(p2; showvalues=[(:iterations, iterations)])
            ProgressMeter.next!(p2; showvalues=[
                (:iterations, iterations),
                (:K, obj_mixture_mcmc[which_first].K[1]),
                (:which_first_temperature, [which_first, temperatures[which_data_to_temp[which_first]]]),
                #(:clusters, countmap(obj_mixture_mcmc[which_first].cluster)),
                (:cc, sort(collect(countmap(obj_mixture_mcmc[which_first].cluster)); by=first))])

            for it in 1:size(temperatures, 1)
                #println(it)
                sampling_mcmc(iterations, obj_graph_mcmc[it], obj_graph_prop[it], obj_mixture_mcmc[it], obj_mixture_prop[it], obj_data_mcmc[it], obj_data_prop[it], obj_prior, temperatures[which_data_to_temp[it]], adaptive_mcmc[it], obj_cohesion_mcmc[it], obj_cohesion_prop[it])
            end


            ### parallel tempering
            if size(temperatures, 1) > 1


                parallel_tempering(iterations, obj_graph_mcmc, obj_graph_prop, obj_mixture_mcmc, obj_mixture_prop, obj_data_mcmc, obj_data_prop, obj_prior, 20, temperatures, which_data_to_temp)

            end

            which_first = findall(x -> x == 1, which_data_to_temp)[1]
            #println("which_first = ", which_first)
        end
        thin_burnin = thin
        ## sample_to_save


        zeta_out[:, iMCMC] = obj_mixture_mcmc[which_first].cluster
        spanning_tree_out[:, :, iMCMC] = obj_graph_mcmc[which_first].graph_st
        w_index_out[:, iMCMC] = (obj_graph_mcmc[which_first].weight_mat[index_node_weight_mat] .- minimum(obj_graph_mcmc[which_first].weight_mat[index_node_weight_mat])) ./ (maximum(obj_graph_mcmc[which_first].weight_mat[index_node_weight_mat]) - minimum(obj_graph_mcmc[which_first].weight_mat[index_node_weight_mat]))



        sep_clust_out[:, :, iMCMC] = obj_mixture_mcmc[which_first].sep_clust
        miss_edge_out[:, :, iMCMC] = obj_mixture_mcmc[which_first].miss_edge
        n_cluster_out[iMCMC] = obj_mixture_mcmc[which_first].K[1]

        mu_out[:, iMCMC] = obj_data_mcmc[which_first].mu
        tau2_out[:, iMCMC] = obj_data_mcmc[which_first].tau2
        sigma2_out[:, iMCMC] = obj_data_mcmc[which_first].sigma2
        rho_out[:, iMCMC] = obj_data_mcmc[which_first].rho
        gp_out[:, :, iMCMC] = obj_data_mcmc[which_first].gp

        predictors_out[:, iMCMC] = obj_graph_mcmc[which_first].predictors
        #K_out = zeros(Int64, sample_to_save)
        #miss_edge_out = zeros(Int64, sample_to_save, 2 * (Kmax - 1))
    end

    return (
        iter_mcmc=iter_mcmc,
        index_node_weight_mat=index_node_weight_mat,
        zeta_out=zeta_out,
        spanning_tree_out=spanning_tree_out,
        w_index_out=w_index_out,
        sep_clust_out=sep_clust_out,
        miss_edge_out=miss_edge_out,
        n_cluster_out=n_cluster_out,
        mu_out=mu_out,
        tau2_out=tau2_out,
        sigma2_out=sigma2_out,
        rho_out=rho_out,
        gp_out=gp_out,
        predictors_out=predictors_out
    )


end



