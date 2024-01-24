%%%%%Posterior permutaion
[A_clu, mu_clu,sigmasq_clu] = post_permute(emit_chain_alex_clus_uninfo,trans_chain_alex_clus_uninfo);
[A_unif,mu_unif,sigmasq_unif] = post_permute(emit_chain_unif,trans_chain_unif);
%%%%%Check predictive likelihood: check on  selected held-out data
z_test = z((T/2+1):T);
indices = randsample(find(z_test==1),100);
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
figure(1)
subplot(1,3,1)
hold on
plot(breakpts, mean(expected_log_density_nonclu,1),'-r','Marker','x')
plot(breakpts, mean(expected_log_density_clu,1),'-b','Marker','x')
ax = gca;
ax.FontSize = 20; 
title("Log Predictive Density (State 1)")
xlabel("Iter Number")
legend("Uniform","TASS")
hold off
%%%%
z_test = z((T/2+1):T);
indices = randsample(find(z_test==2),100);
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
subplot(1,3,2)
hold on
plot(breakpts, mean(expected_log_density_nonclu,1),'-r','Marker','x')
plot(breakpts, mean(expected_log_density_clu,1),'-b','Marker','x')
ax = gca;
ax.FontSize = 20; 
title("Log Predictive Density (State 2)")
xlabel("Iter Number")
legend("Uniform","TASS")
hold off
%%%
%%%%
z_test = z((T/2+1):T);
indices = randsample(find(z_test==3),100);
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
subplot(1,3,3)
hold on
plot(breakpts, mean(expected_log_density_nonclu,1),'-r','Marker','x')
plot(breakpts, mean(expected_log_density_clu,1),'-b','Marker','x')
ax = gca;
ax.FontSize = 20; 
title("Log Predictive Density (State 3)")
xlabel("Iter Number")
legend("Uniform","TASS")
hold off
