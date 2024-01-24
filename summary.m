%%%%Plot logy
plot(y)
yline(-4,'r','X')
yline(-5,'r','M')
yline(-6,'r','C')
title('X-Ray Flux')
disp([mean(y>-4),mean(y>-5),mean(y>-6)]);

%%%%%Posterior permutaion
[A_clu, mu_clu,sigmasq_clu] = post_permute(emit_chain_alex_clus,trans_chain_alex_clus);
[A_unif,mu_unif,sigmasq_unif] = post_permute(emit_chain_unif,trans_chain_unif);
hold on
plot(mu_clu(:,3),'b')
plot(mu_unif(:,3),'r')
title('Traceplot of \mu_3')
hold off

%%%Predictive density conditional on z=4. Treat first 500 as burn-ins.
n_burnin = 500;
pred_y4_clu = zeros(n_mcmc-n_burnin,1);
pred_y4_unif = zeros(n_mcmc-n_burnin,1);
for i = 1 : (n_mcmc-n_burnin)
    pred_y4_clu(i) = normrnd(mu_clu(i,4),sqrt(sigmasq_clu(i,4)));
    pred_y4_unif(i) = normrnd(mu_unif(i,4),sqrt(sigmasq_unif(i,4)));
end
%%%%%conditional on z==4, prob(y_pred>-5|z=4) %%%M-class flare
disp([mean(pred_y4_clu>-5),mean(pred_y4_unif>-5)])
disp([mean(pred_y4_clu>-4),mean(pred_y4_unif>-4)])