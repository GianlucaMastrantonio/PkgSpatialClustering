module PkgSpatClust

using ProgressMeter

using Distributions
using Random
using LinearAlgebra
using PDMats
using StatsBase
using Graphs
using ToggleableAsserts
using SimpleWeightedGraphs
using DelaunayTriangulation
using MetaGraphs
using Clustering
using Combinatorics


#include(joinpath("assert.jl"))
include(joinpath("graphs_types.jl"))
include(joinpath("mixture_types.jl"))
include(joinpath("obs_distributions.jl"))
include(joinpath("obs_types.jl"))

include(joinpath("mcmc.jl"))

include(joinpath("sampling_st.jl"))
include(joinpath("sampling_w.jl"))
include(joinpath("sampling_separators.jl"))
include(joinpath("sampling_ncluster.jl"))
include(joinpath("sampling_mu.jl"))
include(joinpath("sampling_tau2.jl"))
include(joinpath("sampling_sigma2.jl"))
include(joinpath("sampling_rho_and_gp.jl"))
include(joinpath("sampling_cov_par.jl"))

include(joinpath("from_marg_to_full.jl"))


#include(joinpath("sampling_gp.jl"))

export GraphCluter_Vers5, TestMixture_V5, spatial_cluster_v1, GpData_Vers9, PriorsMod1_V3, GpDataMarginalized_Vers9

end # module PkgSpatClust
