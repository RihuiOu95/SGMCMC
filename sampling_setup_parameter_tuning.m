y = detrend_logz'; 
rng(1117676);
t_start = mod(length(y),15)+1;
y = y(t_start:end);
y_test = y((3*1e5+1):end);
y_train = y(1:(3*1e5));
n_latent = 4;
nu = 9;
pi0_true = rand(1,n_latent); 
pi0_true = pi0_true/sum(pi0_true);
emit_dist = string("gaussian");
%%%
L = 2;
B = 5;
eps_clu = 1e-6;
eps_unif = 1e-6;
n_mcmc = 5000;
sgmcmc_param_clu = sgmcmc_parameter(L, B, eps_clu, n_mcmc);
sgmcmc_param_unif = sgmcmc_parameter(L, B, eps_unif, n_mcmc);
N = length(y)/(2*L+1);
%%%%Initial and prior 
trans_param = random_transition_parameter(n_latent, pi0_true);
%emit_param = student_t_emission_parameter([-1,0,1], ones(1,n_latent), nu, emit_dist);
emit_param = gaussian_emission_parameter([-1,0,1], ones(1,n_latent), emit_dist);
mb_size = 10;
mu_prior = gaussian_prior(zeros(1,n_latent), 10*ones(1,n_latent));
sigmasq_prior = IG_prior(3*ones(1,n_latent), 10*ones(1,n_latent));
emit_prior = gaussian_emission_prior(mu_prior, sigmasq_prior);
%%%%%Weights
clear wt_alex_clus wt_unif
eta = 1e-5;
[wt_alex_clus, center, xisum] = weights_alex_clus_leaveoneout(y, n_latent, N, L, eta);
wt_unif = weights_unif(n_latent, N);
A_prior = reshape(xisum, n_latent, n_latent);
prior_prob_mat = reshape(xisum,n_latent,n_latent);
prior_prob_mat = prior_prob_mat ./ sum(prior_prob_mat,1);
%trans_prior_info = dirichlet_prior(1e7*prior_prob_mat);
trans_prior_uninfo = dirichlet_prior(1*prior_prob_mat);
%%%%%%%
trans_param_init = random_transition_parameter(n_latent, pi0_true);
%emit_param_init = student_t_emission_parameter(zeros([1,n_latent]), 0.1*ones(1,n_latent),nu, emit_dist);
emit_param_init = gaussian_emission_parameter(zeros([1,n_latent]), 0.1*ones(1,n_latent), emit_dist);
clear trans_chain_alex_clus emit_chain_alex_clus trans_chain_alex_orac emit_chain_alex_orac 
clear trans_chain_unif emit_chain_unif
%[trans_chain_alex_clus_info, emit_chain_alex_clus_info, A_hat_chain_clu, A_chain_clu] = CSG_MCMC_IS(y, wt_alex_clus, mb_size, emit_param_init, ... 
                                                           %trans_param_init, sgmcmc_param_clu, emit_prior, trans_prior_info);
[trans_chain_alex_clus_uninfo, emit_chain_alex_clus_uninfo, A_hat_chain_clu, A_chain_clu] = CSG_MCMC_IS(y, wt_alex_clus, mb_size, emit_param_init, ... 
                                                           trans_param_init, sgmcmc_param_clu,1e-6, emit_prior, trans_prior_uninfo);
[trans_chain_unif, emit_chain_unif, A_hat_chain_unif] = CSG_MCMC_IS(y, wt_unif, mb_size, emit_param_init, ... 
                                                      trans_param_init, sgmcmc_param_unif,1e-6, emit_prior, trans_prior_uninfo);