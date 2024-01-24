A_true_clu = reshape(trans_chain_alex_clus_uninfo.A(n_mcmc,:,:),n_latent,n_latent);
A_true_nonclu = reshape(trans_chain_unif.A(n_mcmc,:,:),n_latent,n_latent);
%%%%
%emit_param_true_clu = student_t_emission_parameter(emit_chain_alex_clus.mu(n_mcmc,:),emit_chain_alex_clus.sigmasq(n_mcmc,:),nu, emit_dist);
%emit_param_true_nonclu = student_t_emission_parameter(emit_chain_unif.mu(n_mcmc,:),emit_chain_unif.sigmasq(n_mcmc,:),nu, emit_dist);
emit_param_true_clu = gaussian_emission_parameter(emit_chain_alex_clus_uninfo.mu(n_mcmc,:),emit_chain_alex_clus_uninfo.sigmasq(n_mcmc,:), emit_dist);
emit_param_true_nonclu = gaussian_emission_parameter(emit_chain_unif.mu(n_mcmc,:),emit_chain_unif.sigmasq(n_mcmc,:), emit_dist);
%%%%
[z_pred_clu, y_pred_clu] = gen_rep([1,1], length(y), 1,  A_true_clu, emit_param_true_clu);
[z_pred_nonclu, y_pred_nonclu] = gen_rep([1,1], length(y), 1,  A_true_nonclu, emit_param_true_nonclu);
%%%%
figure
subplot(3,1,1)
plot(y_pred_clu)
title("Simulated Series---Clustering")
ylim([-1.5 2.5])
subplot(3,1,2)
plot(y_pred_nonclu)
ylim([-1.5 2.5])
title("Simulated Series---Uniform")
subplot(3,1,3)
plot(y)
ylim([-1.5 2.5])
title("Data")
%%%%%Check predictive likelihood: check on  selected held-out data
cutoff = 1.0;
indices = identify_brust_region(y_test,cutoff,100);
breakpts = (100:100:3000);
expected_log_density_clu = zeros([length(indices),length(breakpts)]);
expected_log_density_nonclu = zeros([length(indices), length(breakpts)]);
n_test = 10;
for i = 1 : length(indices)
    t_start = indices(i)-1;
    y_held = y_test(indices(i):(indices(i)+n_test));
    buffer = 10;
    expected_log_density_nonclu(i,:) = lppd(y_held, y_test, trans_chain_unif, emit_chain_unif,  t_start, emit_dist, breakpts, buffer, nu);
    expected_log_density_clu(i,:) = lppd(y_held, y_test, trans_chain_alex_clus_uninfo, emit_chain_alex_clus_uninfo,  t_start, emit_dist, breakpts, buffer, nu);
end

figure
%subplot(2,1,1)
hold on
plot(breakpts, mean(expected_log_density_nonclu,1),'-r','Marker','x')
plot(breakpts, mean(expected_log_density_clu,1),'-b','Marker','x')
ax = gca;
ax.FontSize = 16; 
%title("Holding Out Rare Events-Log Predictive Density")
xlabel("Iter Number")
ylabel("Log Predictive Density")
legend("Uniform","TASS")
hold off
%%%%% Check Coverage
n_burnin = 2000;
n_step_ahead = 5;
y_pred_samps_nonclu = zeros([length(indices),n_mcmc-n_burnin, n_step_ahead]);
y_pred_samps_clu = zeros([length(indices),n_mcmc-n_burnin, n_step_ahead]);
for i = 1 : length(indices)
    t_start = indices(i)-1;
    buffer = 5;
    y_pred_samps_nonclu(i,:,:) = gen_pred_obs(y, trans_chain_unif, emit_chain_unif,  t_start, emit_dist, n_burnin, n_step_ahead, buffer, nu);
    y_pred_samps_clu(i,:,:) = gen_pred_obs(y, trans_chain_alex_clus_uninfo, emit_chain_alex_clus_uninfo,  t_start, emit_dist, n_burnin, n_step_ahead, buffer, nu);
end

confin_lvls = 0.8 : 0.03 : 0.95;
n_lvl = length(confin_lvls);
cover_prob_clu = zeros([n_lvl,1]);
cover_prob_nonclu = zeros([n_lvl,1]);
for j = 1:n_lvl
    confin_lvl = confin_lvls(j);
    lower = (1-confin_lvl)/2;
    upper = (1+confin_lvl)/2;
    clu_cover_counts = 0;
    nonclu_cover_counts = 0;
    for i = 1 : length(indices)
        t_start = indices(i)-1;
        y_held = y(indices(i):(indices(i)+n_step_ahead-1))';
        confint_clu = quantile(reshape(y_pred_samps_clu(i,:,:),[n_mcmc-n_burnin, n_step_ahead]),[upper,lower]);
        confint_nonclu = quantile(reshape(y_pred_samps_nonclu(i,:,:),[n_mcmc-n_burnin, n_step_ahead]),[upper,lower]);
        clu_cover_counts = clu_cover_counts+ sum((y_held<confint_clu(1,:)) & (y_held>confint_clu(2,:)));
        nonclu_cover_counts = nonclu_cover_counts+ sum((y_held<confint_nonclu(1,:)) & (y_held>confint_nonclu(2,:)));
    end
    cover_prob_clu(j) = clu_cover_counts / n_step_ahead / length(indices);
    cover_prob_nonclu(j) = nonclu_cover_counts / n_step_ahead / length(indices);
end


figure
hold on
plot(confin_lvls, confin_lvls, '--magenta')
plot(confin_lvls, cover_prob_nonclu,'-r','Marker','x')
plot(confin_lvls, cover_prob_clu,'-b','Marker','x')
ax = gca;
ax.FontSize = 16; 
%title("Holding Out Rare Events-Log Predictive Density")
xlabel("Theroatical Coverage Probability")
ylabel("Actual Coverage Probability")
legend("Theoretical", "Uniform","Clustering")
hold off



%%%%%
%%%%%%%%
%%%%%Step 2: Check on random held-out test
indices = randsample(100:length(z)-100,60);
breakpts = (100:200:2000);
expected_log_density_clu_unifheld = zeros([length(indices),length(breakpts)]);
expected_log_density_nonclu_unifheld = zeros([length(indices), length(breakpts)]);
n_test = 10;
for i = 1 : length(indices)
    t_start = indices(i)-1;
    y_held = y(indices(i):(indices(i)+n_test));
    buffer = 10;
    expected_log_density_nonclu_unifheld(i,:) = lppd(y_held, y, trans_chain_unif, emit_chain_unif,  t_start, emit_dist, breakpts, buffer);
    expected_log_density_clu_unifheld(i,:) = lppd(y_held, y, trans_chain_alex_clus_uninfo, emit_chain_alex_clus_uninfo,  t_start, emit_dist, breakpts, buffer);
end

figure
%subplot(2,1,1)
hold on
scatter(breakpts, mean(expected_log_density_nonclu,1),'r','Marker','x')
scatter(breakpts, mean(expected_log_density_clu,1),'b')
title("Holding Out Rare Events-Log Predictive Density")
xlabel("Iter Number")
legend("Uniform","Clutering")
hold off

subplot(2,1,2)
hold on
scatter(breakpts, mean(expected_log_density_nonclu_unifheld,1),'r')
scatter(breakpts, mean(expected_log_density_clu_unifheld,1),'b')
title("Holding Out Rare Events")
xlabel("Iter Number")
legend("Uniform","Clutering")
hold off

