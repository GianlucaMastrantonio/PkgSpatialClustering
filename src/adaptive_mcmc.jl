abstract type AdaptiveMCMC end

struct MultiType1_V3 <: AdaptiveMCMC

  dimension::Int64
  par_acc::Vector{Float64}
  par_prop::Vector{Float64}
  mean_vec::Vector{Float64}
  cov_mat::Matrix{Float64}
  lambda::Vector{Float64}

  mean_init::Vector{Float64}
  cov_init::Matrix{Float64}
  lambda_init::Vector{Float64}

  diag_mat::Matrix{Float64}

  target_alpha::Float64
  accept_rate::Vector{Float64}

  iter_start::Int64
  iter_end::Int64

  par_a::Float64
  par_b::Float64

  number_iter_adaptive::Vector{Int64}

  function MultiType1_V3(
    dimension::Int64,
    par_acc::Vector{Float64},
    par_prop::Vector{Float64},
    mean_vec::Vector{Float64},
    cov_mat::Matrix{Float64},
    lambda::Vector{Float64},

    mean_init::Vector{Float64},
    cov_init::Matrix{Float64},
    lambda_init::Vector{Float64},

    diag_mat::Matrix{Float64},
    target_alpha::Float64,
    accept_rate::Vector{Float64},
    iter_start::Int64,
  iter_end::Int64,
  par_a::Float64,
  par_b::Float64,
    number_iter_adaptive::Vector{Int64}
  )

    new(dimension, par_acc, par_prop, mean_vec, cov_mat, lambda, mean_init, cov_init, lambda_init, diag_mat, target_alpha, accept_rate, iter_start, iter_end, par_a, par_b, number_iter_adaptive)

  end

end

function MultiType1_V3(; parameters::Vector{Float64},
  cov_init::Matrix{Float64},
  lambda_init::Vector{Float64},
  diag_mat::Matrix{Float64},
  target_alpha::Float64,
  iter_start::Int64, 
  iter_end::Int64,
  par_a::Float64,
  par_b::Float64)

  dimension::Int64 = size(parameters, 1)
  mean_init::Vector{Float64} = [rand(Normal(parameters[i], diag_mat[i,i])) for i in 1:dimension]
  mean_vec::Vector{Float64} = deepcopy(mean_init)
  cov_mat::Matrix{Float64} = deepcopy(cov_init)
  lambda::Vector{Float64} = deepcopy(lambda_init)
  accept_rate::Vector{Float64} = [0.0]

  par_acc::Vector{Float64} = deepcopy(parameters)
  par_prop::Vector{Float64} = deepcopy(parameters)

  MultiType1_V3(dimension, par_acc, par_prop, mean_vec, cov_mat, lambda, mean_init, deepcopy(cov_init), deepcopy(lambda_init), deepcopy(diag_mat), target_alpha, accept_rate, iter_start, iter_end, par_a, par_b, [0])


end


function sample_amcmc(amcmc::MultiType1_V3)

  amcmc.par_prop .= amcmc.lambda[1]^0.5 .* cholesky(amcmc.cov_mat + amcmc.diag_mat).L * rand(Normal(), amcmc.dimension) .+ amcmc.par_acc

end

function update_amcmc(amcmc::MultiType1_V3, iteration::Int64, accept::Bool, alpha_mh::Float64)


  if iteration >= amcmc.iter_start & iteration <= amcmc.iter_end
    amcmc.number_iter_adaptive[1] += 1
    if accept == true
      amcmc.par_acc .= amcmc.par_prop
      amcmc.accept_rate[1] = (amcmc.accept_rate[1] * (amcmc.number_iter_adaptive[1] - 1) + 1.0) / amcmc.number_iter_adaptive[1]
    else
      amcmc.accept_rate[1] = (amcmc.accept_rate[1] * (amcmc.number_iter_adaptive[1] - 1) + 0.0) / amcmc.number_iter_adaptive[1]
    end
    gamma::Float64 = amcmc.par_a / (amcmc.par_b + iteration) 
    amcmc.lambda[1] = exp(log(amcmc.lambda[1]) + gamma* (alpha_mh - amcmc.target_alpha))
    amcmc.mean_vec .+= gamma * (amcmc.par_acc - amcmc.mean_vec)
    amcmc.cov_mat .+= gamma * ((amcmc.par_acc - amcmc.mean_vec) * transpose(amcmc.par_acc.- amcmc.mean_vec) - amcmc.cov_mat)

    #if is_pos_def(amcmc.cov_mat)


    #end
    

  end

  


end