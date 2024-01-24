function weights = weights_true(z, n_latent, K, L, rare_state, eta) 
    clear w weights_trans wt_emit weights_emit
    ll = 2*L+1;
    for k = 1:K
        idx = z(((k-1)*ll+1):(k*ll));
        w(k) = sum(idx==rare_state) + eta;
    end
    
    for j = 1:n_latent 
        for k = 1:n_latent 
            clear wt
            if j == rare_state || k == rare_state
                wt = w;
            else 
                wt = ones(1,K);
            end
            weights_trans{j,k} = wt/sum(wt);
        end
    end
    for j = 1:n_latent 
        clear wt
        if j == rare_state
            wt = w;
        else 
            wt = ones(1,K);
        end
        wt_emit{j} = wt/sum(wt);
    end
    
    weights_emit = weight_emit_gaussian(wt_emit, wt_emit);

    weights = weight(weights_trans, weights_emit);
end