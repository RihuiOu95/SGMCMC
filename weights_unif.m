function weights = weights_unif(n_latent, K) 
    wt = ones(1,K);
    clear weights_trans wt_emit weights_emit
    for j = 1:n_latent 
        for k = 1:n_latent 
            weights_trans{j,k} = wt/sum(wt);
        end
    end
    for j = 1:n_latent 
        wt_emit{j} = wt/sum(wt);
    end
    weights_emit = weight_emit_gaussian(wt_emit, wt_emit);
    weights = weight(weights_trans, weights_emit);
end