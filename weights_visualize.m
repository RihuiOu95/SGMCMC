rng(11176767);
A = [.99 .01;.01 .99];
mu = [-20,20];
sigmasq = [1,1];
pi0 = [.5, .5];
trans_param = transition_parameter(pi0, A, A);
emit_param = gaussian_emission_parameter(mu, sigmasq, "gaussian");
T = 1e4;
[z,y] = simHMM(trans_param, emit_param, T);
y_test = y((T/2+1):T);
y = y(1:T/2);
%%%%%%
%t_start = mod(length(y),15)+1;
%y = y(t_start:end);
n_latent = 3;
pi0_true = rand(1,n_latent); 
pi0_true = pi0_true/sum(pi0_true);
emit_dist = string("gaussian");
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
[wt_clus, center, xisum] = weights_alex_clus_leaveoneout(y, n_latent, N, L, eta);
wt_unif = weights_unif(n_latent, N);
%%%%%Plot the weights
figure
hold on
plot(log10(wt_clus.emit.mu{1,3}))
title("Log Predictive Weights")
xlabel("T")
ax = gca;
ax.FontSize = 16; 
saveas(gcf,'results/weights_nonrare.png')
hold off