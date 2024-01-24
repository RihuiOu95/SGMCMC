function [clus_res,unif_res] = sampling_setup_gradient_mse(T, L, B, std)
    rng(1117676);
    A = [.99 .005 .4950;.005 .99 .4950;.005 .005 .0100];
    mu = [-20, 0, 20];
    sigmasq = [1,1,1];
    pi0 = [1/3 1/3 1/3];
    trans_param = transition_parameter(pi0, A, A);
    emit_param = gaussian_emission_parameter(mu, sigmasq, "gaussian");
    %%%
    T = T * 2;
    [z,y] = simHMM(trans_param, emit_param, T);
    y_test = y((T/2+1):T);
    y = y(1:T/2);
    n_latent = 3;
    pi0_true = rand(1,n_latent); 
    pi0_true = pi0_true/sum(pi0_true);
    emit_dist = string("gaussian");
    %L = 7;
    %B = 5;
    eps_clu = 1e-6;
    eps_unif = 1e-6;
    n_mcmc = 5000;
    sgmcmc_param_clu = sgmcmc_parameter(L, B, eps_clu, n_mcmc);
    sgmcmc_param_unif = sgmcmc_parameter(L, B, eps_unif, n_mcmc);
    N = length(y)/(2*L+1);
    %%%%Initial and prior 
    trans_param = random_transition_parameter(n_latent, pi0_true);
    emit_param = gaussian_emission_parameter([-20, 0, 20] + std, ones(1,n_latent), emit_dist);
    mb_size = 1;
    mu_prior = gaussian_prior(zeros(1,n_latent), 10*ones(1,n_latent));
    sigmasq_prior = IG_prior(3*ones(1,n_latent), 10*ones(1,n_latent));
    emit_prior = gaussian_emission_prior(mu_prior, sigmasq_prior);
    %%%%%Weights
    clear wt_alex_clus wt_unif
    eta = 1e-5;
    [wt_alex_clus, center, xisum] = weights_alex_clus_leaveoneout(y, n_latent, N, L, eta);
    wt_unif = weights_alex_clus(y, n_latent, N, L, eta);
    A_prior = reshape(xisum, n_latent, n_latent);
    prior_prob_mat = reshape(xisum,n_latent,n_latent);
    prior_prob_mat = prior_prob_mat ./ sum(prior_prob_mat,1);
    trans_prior_uninfo = dirichlet_prior(zeros([n_latent,n_latent])+1);
    trans_param_init = random_transition_parameter(n_latent, pi0_true);
    emit_param_init = gaussian_emission_parameter(zeros([1,n_latent]), 5*ones(1,n_latent), emit_dist);
    clear trans_chain_alex_clus emit_chain_alex_clus trans_chain_alex_orac emit_chain_alex_orac 
    clear trans_chain_unif emit_chain_unif
    %emit_param.eval_full_grad(y, trans_param, sgmcmc_param_clu, emit_prior)
    %emit_param.eval_stoch_grad_IS(y, wt_alex_clus,  mb_size, trans_param, sgmcmc_param_clu, emit_prior)
    %emit_param.eval_stoch_grad_IS(y, wt_unif,  mb_size, trans_param, sgmcmc_param_clu, emit_prior)

    n_rep = 100;
    trans_param_full = trans_param.eval_grad_full(y, emit_param, sgmcmc_param_clu, trans_prior_uninfo);
    
    % find the rare label;
    rare_label = find(abs(center - 20) == min(abs(center - 20)));
    ground_truth = trans_param_full.A_grad(rare_label,:);
    unif_grad_error = zeros(n_rep, n_latent);
    clus_grad_error = zeros(n_rep, n_latent);
    for i = 1 : n_rep
        emit_param_clus = trans_param.eval_stoch_grad_IS(y, wt_alex_clus,  mb_size, emit_param, sgmcmc_param_clu, trans_prior_uninfo);
        emit_param_unif = trans_param.eval_stoch_grad_IS(y, wt_unif,  mb_size, emit_param, sgmcmc_param_clu, trans_prior_uninfo);
        clus_grad_error(i, :) = emit_param_clus.A_grad(rare_label, :);
        unif_grad_error(i, :) = emit_param_unif.A_grad(rare_label, :); 
    end
    
    clus_grad_error = clus_grad_error - ground_truth;
    unif_grad_error = unif_grad_error - ground_truth;
    clus_res = sqrt(sum(sum(clus_grad_error.^2)) / n_rep);
    unif_res = sqrt(sum(sum(unif_grad_error.^2)) / n_rep);
end

