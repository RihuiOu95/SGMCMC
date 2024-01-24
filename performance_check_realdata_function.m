function [unif_log_pred, clu_log_pred] = tuning_eval(trans_chain_alex_clus_uninfo, trans_chain_unif, emit_chain_alex_clus_uninfo, emit_chain_unif)
    cutoff = 0.8;
    indices = identify_brust_region(y_test,cutoff,10);
    breakpts = n_mcmc;
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
    unif_log_pred = mean(expected_log_density_nonclu);
    clu_log_pred = mean(expected_log_density_clu);
end


