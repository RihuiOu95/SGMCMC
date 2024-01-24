%%%%%Posterior permutaion
[A_clu, mu_clu,sigmasq_clu] = post_permute(emit_chain_alex_clus_uninfo,trans_chain_alex_clus_uninfo);
[A_unif,mu_unif,sigmasq_unif] = post_permute(emit_chain_unif,trans_chain_unif);

%%%%Check mu_3 - mu_3^{true}
mu_true = 20;
mu_clu_error = abs(mu_clu(:,3)-mu_true);
mu_unif_error = abs(mu_unif(:,3)-mu_true);
figure(1)
%title('|\mu_3 - \mu_3^{true}|')
plot(mu_clu_error,'b')
hold on 
plot(mu_unif_error,'r')
xlim([0, n_mcmc])
ax = gca;
ax.FontSize = 16; 
%title("Holding Out Rare Events-Log Predictive Density")
xlabel("Iter Number")
ylabel("Error")
legend("TASS","Uniform")
hold off

%%%%Check sigma^2_3 - sigma_3^{true}
sigmasq_true = 1;
sigmasq_clu_error = abs(sigmasq_clu(:,3)-sigmasq_true);
sigmasq_unif_error = abs(sigmasq_unif(:,3)-sigmasq_true);
figure(1)
plot(sigmasq_clu_error,'b')
hold on 
plot(sigmasq_unif_error,'r')
xlim([0, n_mcmc])
ax = gca;
ax.FontSize = 16; 
%title("Holding Out Rare Events-Log Predictive Density")
xlabel("Iter Number")
ylabel("Error")
legend("TASS","Uniform")
hold off







%%%%%Check predictive likelihood: check on  selected held-out data
z_test = z((T/2+1):T);
indices = randsample(find(z_test==3),200);
breakpts = (100:200:5000);
expected_log_density_clu = zeros([length(indices),length(breakpts)]);
expected_log_density_nonclu = zeros([length(indices), length(breakpts)]);
n_test = 0;
nu = 9;
for i = 1 : length(indices)
    t_start = indices(i)-1;
    y_held = y_test(indices(i):(indices(i)+n_test));
    buffer = 1;
    expected_log_density_nonclu(i,:) = lppd(y_held, y_test, trans_chain_unif, emit_chain_unif,  t_start, emit_dist, breakpts, buffer, nu);
    expected_log_density_clu(i,:) = lppd(y_held, y_test, trans_chain_alex_clus_uninfo, emit_chain_alex_clus_uninfo,  t_start, emit_dist, breakpts, buffer, nu);
end
%%%%Plot
figure
hold on
plot(breakpts, mean(expected_log_density_nonclu,1),'-r','Marker','x')
plot(breakpts, mean(expected_log_density_clu,1),'-b','Marker','x')
ax = gca;
ax.FontSize = 16; 
title("Log Predictive Density (Simulated Data)")
xlabel("Iter Number")
legend("Uniform","TASS")
hold off
