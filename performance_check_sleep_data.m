

%%%%%Check predictive likelihood: check on  selected held-out data
cutoff = 3.8;
indices = identify_brust_region(y_test,cutoff, 100, true);
breakpts = (300:300:n_mcmc);
expected_log_density_clu = zeros([length(indices),length(breakpts)]);
expected_log_density_nonclu = zeros([length(indices), length(breakpts)]);
n_test = 0;
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
title("Log Predictive Density (Sleep Data)")
xlabel("Iter Number")
ylabel("Log Predictive Density")
legend("Uniform","TASS")
hold off
