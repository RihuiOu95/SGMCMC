%This script is mainly to evaluate the predictive performance.
n_test = 10; %length of each held out dataset
n_rep = 20; %total numbers of replicate dataset
%%%%%Generate n_rep replicate data
[z_rep_1,y_rep_1] = gen_rep(z_1, n_test, n_rep,  A_true, emit_param_true);
%%%%%Calculate Log-predictive density
%%%%%First obtain trans_chain_unif_1...... from sampler.m
t_start = T;
lpd_unif_1 = lppd(y_rep_1, y_1,  trans_chain_unif_1, emit_chain_unif_1, sgmcmc_param, t_start, emit_dist,(100:100:1000));
lpd_cluster_1 = lppd(y_rep_1, y_1,  trans_chain_cluster_1, emit_chain_cluster_1, sgmcmc_param, t_start, emit_dist,(100:100:1000));
lpd_true_1 = lppd(y_rep_1, y_1,  trans_chain_true_1, emit_chain_true_1, sgmcmc_param, t_start, emit_dist,(100:100:1000));
%%%%%Plot
title('10-step-log-predictive denity')
hold on
plot((100:100:1000),lpd_unif_1)
plot((100:100:1000),lpd_cluster_1)
plot((100:100:1000),lpd_true_1)
xlabel("MCMC iteration")
ylabel("log-predictive")
legend('uniform','cluster','oracle')
hold off
%%%%%%Relative Performance---10%
[z_rep_10,y_rep_10] = gen_rep(z_10, n_test, n_rep,  A_true_10pct, emit_param_true);
lpd_unif_10 = lppd(y_rep_10, y_10,  trans_chain_unif_10, emit_chain_unif_10, sgmcmc_param, t_start, emit_dist,1000);
lpd_cluster_10 = lppd(y_rep_10, y_10,  trans_chain_cluster_10, emit_chain_cluster_10, sgmcmc_param, t_start, emit_dist,1000);

%%%%%%Relative Performance---balance
[z_rep_bd,y_rep_bd] = gen_rep(z_balance, n_test, n_rep,  A_true_balance, emit_param_true);
lpd_unif_bd = lppd(y_rep_bd, y_balance,  trans_chain_unif_balance, emit_chain_unif_balance, sgmcmc_param, t_start, emit_dist,1000);
lpd_cluster_bd = lppd(y_rep_bd, y_balance,  trans_chain_cluster_balance, emit_chain_cluster_balance, sgmcmc_param, t_start, emit_dist,1000);
%%%%%Part II Calculate ||A-A_true||
%%%%%Permute posterior samples
[A_unif_sorted,mu_unif_sorted,sigmasq_unif_sorted] = post_permute(emit_chain_unif_1,trans_chain_unif_1);
[A_cluster_sorted,mu_cluster_sorted,sigmasq_cluster_sorted] = post_permute(emit_chain_cluster_1,trans_chain_cluster_1);
[A_true_sorted,mu_true_sorted,sigmasq_true_sorted] = post_permute(emit_chain_true_1,trans_chain_true_1);
%%%%%

clear A_error_true A_error_cluster A_error_unif
for i = 1:n_mcmc 
    A_error_true(i) = norm(reshape(A_true_sorted(i,:,:),n_latent,n_latent) - A_true_1pct);
    A_error_cluster(i) = norm(reshape(A_cluster_sorted(i,:,:),n_latent,n_latent) - A_true_1pct);
    A_error_unif(i) = norm(reshape(A_unif_sorted(i,:,:),n_latent,n_latent) - A_true_1pct);
end

hold on
plot(A_error_unif)
title("Uniform sub-sampling")
ylabel("|A-A_{true}|")
xlabel("MCMC iteration")
grid on

plot(A_error_cluster)
title("Importance sampling with clustered y")
ylabel("|A-A_{true}|")
xlabel("MCMC iteration")
grid on
hold off



%%%%%
disp([lpd_unif,lpd_cluster,lpd_true])

%%%%Things to be done immediately:
%1. Adapt the lppd.m to generate a log-pred density vs iterations
%plot;

%2. Write post_permute.m in an OOP way;

%3. Use more performance metrics--hold out some rare states to measure how
%well this method can perdict rare state parameters.(later)

%4. Relative Performance;

%5. Update the draft

%%%%Question
% As rare states becoecome rarer, how to measure relative performance?
% In this case, the uniform subsampling does not converge at all, so it
% can't be worse in some sense---like throwing a dice from 1-3 to estimate
% rare states. However, the clustering subsampling may presumbably be worse
% and thus we'll see decline in relative performance.