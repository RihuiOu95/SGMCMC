%%%%%Check predictive likelihood: check on  selected held-out data
z_test = z((T/2+1):T);
indices = randsample(find(z_test==3), 100);
breakpts = (n_mcmc/10):(n_mcmc/10):n_mcmc;
breakpts_full = (n_mcmc_full/10):(n_mcmc_full/10):n_mcmc_full;
expected_log_density_clu = zeros([length(indices),length(breakpts)]);
expected_log_density_nonclu = zeros([length(indices), length(breakpts)]);
expected_log_density_full = zeros([length(indices), length(breakpts_full)]);
n_test = 0;
nu = 9;
for i = 1 : length(indices)
    t_start = indices(i)-1;
    y_held = y_test(indices(i):(indices(i)+n_test));
    buffer = 5;
    expected_log_density_nonclu(i,:) = lppd(y_held, y_test, trans_chain_unif, emit_chain_unif,  t_start, emit_dist, breakpts, buffer, nu);
    expected_log_density_clu(i,:) = lppd(y_held, y_test, trans_chain_alex_clus_uninfo, emit_chain_alex_clus_uninfo,  t_start, emit_dist, breakpts, buffer, nu);
    expected_log_density_full(i,:) = lppd(y_held, y_test, trans_chain_full, emit_chain_full,  t_start, emit_dist, breakpts_full, buffer, nu);
end
%%%%Plot

figure
hold on
plot(log10(runtimes_full), mean(expected_log_density_full,1),'-k','Marker','x', 'LineWidth',1.5)
plot(log10(runtimes_unif), mean(expected_log_density_nonclu,1),'-r','Marker','x', 'LineWidth',1.5)
plot(log10(runtimes_clu), mean(expected_log_density_clu,1),'-b','Marker','x', 'LineWidth',1.5)
ax = gca;
ax.FontSize = 18; 
set(gca,'XTickLabel',{'10^{1.5}', '10^2', '10^{2.5}', '10^3', '10^{3.5}', '10^4', '10^{4.5}', '10^5'})
title("Log Predictive Density")
xlabel("CPU Times (secs)")
legend("Full","Uniform", "TASS")
hold off