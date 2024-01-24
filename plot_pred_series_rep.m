figure;
for i = 1:3
    iter_num = n_mcmc - (i-1)*500;
    A_true_clu = reshape(trans_chain_alex_clus_uninfo.A(iter_num,:,:),n_latent,n_latent);
    emit_param_true_clu = gaussian_emission_parameter(emit_chain_alex_clus_uninfo.mu(iter_num,:),emit_chain_alex_clus_uninfo.sigmasq(iter_num,:), emit_dist);
    [z_pred_clu, y_pred_clu] = gen_rep([1,1], length(y), 1,  A_true_clu, emit_param_true_clu);
    subplot(3,2,2*i-1)
    plot(y_pred_clu);
    ax = gca;
    ax.FontSize = 22; 
    ylim([-1.5 2.5])
    %rare_prob = cal_stationary(trans_chain_alex_clus_uninfo,iter_num);
    %rare_prob = rare_prob(3);
    %rare_prob = round(rare_prob,4) * 100;
    %title(strcat("Predicted Series of TASS at ",string(iter_num),"th iteration"))     
end

%%%Uniform
for i = 1:3
    iter_num = n_mcmc - (i-1)*500;
    A_true_clu = reshape(trans_chain_unif.A(iter_num,:,:),n_latent,n_latent);
    emit_param_true_nonclu = gaussian_emission_parameter(emit_chain_unif.mu(iter_num,:),emit_chain_unif.sigmasq(iter_num,:), emit_dist);
    [z_pred_clu, y_pred_clu] = gen_rep([1,1], length(y), 1,  A_true_clu, emit_param_true_nonclu);
    subplot(3,2,2*i)
    plot(y_pred_clu)
    ylim([-1.5 2.5])
    rare_prob = cal_stationary(trans_chain_unif,iter_num);
    rare_prob = rare_prob(3);
    rare_prob = round(rare_prob,4) * 100;
    %title(strcat("Predicted Series of Uniform Subsampler at ",string(iter_num),"th iteration"))
    ax = gca;
    ax.FontSize = 22;  
end

%%%%
subplot(8,2,15);
plot(y)
ylim([-1.5 2.5])
title("Data")