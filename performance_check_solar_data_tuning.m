Ls = [2, 7];
Bs = [5,10];
mb_sizes = [2, 5, 10];
n_mcmc = 200;
res = zeros(length(Ls),length(Bs),length(mb_sizes),2);
%%%%%%%%%%
n_latent = 4;
pi0_true = rand(1,n_latent); 
pi0_true = pi0_true/sum(pi0_true);
emit_dist = string("gaussian");
%%%%%set up prior
%%%%Initial and prior 
trans_param = random_transition_parameter(n_latent, pi0_true);
emit_param = gaussian_emission_parameter([-1,0,1], ones(1,n_latent), emit_dist);
mu_prior = gaussian_prior(zeros(1,n_latent), 10*ones(1,n_latent));
sigmasq_prior = IG_prior(3*ones(1,n_latent), 10*ones(1,n_latent));
emit_prior = gaussian_emission_prior(mu_prior, sigmasq_prior);
%%%%%%%
trans_param_init = random_transition_parameter(n_latent, pi0_true);
emit_param_init = gaussian_emission_parameter(zeros([1,n_latent]), 0.1*ones(1,n_latent), emit_dist);
%%%%
%%%%%test tuning parameters%%%%%
for i = 1 : length(Ls)
    L = Ls(i);
    %%%%%Weights for different Ls
    clear wt_alex_clus wt_unif
    eta = 1e-5;
    N = length(y)/(2*L+1);
    [wt_alex_clus, center, xisum] = weights_alex_clus_leaveoneout(y, n_latent, N, L, eta);
    wt_unif = weights_unif(n_latent, N);
    A_prior = reshape(xisum, n_latent, n_latent);
    prior_prob_mat = reshape(xisum,n_latent,n_latent);
    prior_prob_mat = prior_prob_mat ./ sum(prior_prob_mat,1);
    trans_prior_uninfo = dirichlet_prior(1*prior_prob_mat);
    %%%%%%%
    for j = 1 : length(Bs)
        B = Bs(j);
        for k = 1 : length(mb_sizes)
            %%%%%%
            mb_size = mb_sizes(k);
            %%%%start_doing_MCMC
            disp('Doing MCMC on tuning paramters of')
            disp([L, B, mb_size])
            eps_clu = 1e-7;
            eps_unif = 1e-7;
            sgmcmc_param_clu = sgmcmc_parameter(L, B, eps_clu, n_mcmc);
            sgmcmc_param_unif = sgmcmc_parameter(L, B, eps_unif, n_mcmc);
            clear trans_chain_alex_clus emit_chain_alex_clus trans_chain_alex_orac emit_chain_alex_orac 
            clear trans_chain_unif emit_chain_unif
            [trans_chain_alex_clus_uninfo, emit_chain_alex_clus_uninfo, A_hat_chain_clu, A_chain_clu] = CSG_MCMC_IS(y, wt_alex_clus, mb_size, emit_param_init, ... 
                                                       trans_param_init, sgmcmc_param_clu,1e-6, emit_prior, trans_prior_uninfo);
            [trans_chain_unif, emit_chain_unif, A_hat_chain_unif] = CSG_MCMC_IS(y, wt_unif, mb_size, emit_param_init, ... 
                                                  trans_param_init, sgmcmc_param_unif,1e-6, emit_prior, trans_prior_uninfo);
            [res(i, j, k, 1), res(i, j, k, 2)] = tuning_eval(trans_chain_alex_clus_uninfo, trans_chain_unif, emit_chain_alex_clus_uninfo, emit_chain_unif, y_test, n_mcmc);
        end
    end
end