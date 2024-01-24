function lpd = log_pred_dens(y, trans_chain, emit_chain, sgmcmc_param, t_start, chunk_len, emit_dist) 
    n_mcmc = size(trans_chain.A,1);
    n_latent = size(trans_chain.A,2);
    lpd = zeros(1,n_mcmc);
    for s = 1:n_mcmc 
        A = reshape(trans_chain.A(s,:,:), n_latent,n_latent);
        trans_param = transition_parameter(ones(1,n_latent)/n_latent, A, A);
        emit_param = gaussian_emission_parameter(emit_chain.mu(s,:), emit_chain.sigmasq(s,:), emit_dist);
        [pipred,qpred] = approx_pred(y, t_start, sgmcmc_param, trans_param, emit_param);
        pred = ones(1,n_latent);
        for t = (t_start+chunk_len):-1:(t_start+1)
            pred = pred*emit_param.calP(y(t))*trans_param.A;
        end
        lpd(s) = log(pred*qpred);
    end
end