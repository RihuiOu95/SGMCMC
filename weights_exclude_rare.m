function weights = weights_exclude_rare(z, n_latent, K, L, rare_state) 
    clear w weights_trans weight_emit 
    ll = 2*L+1;
    for k = 1:K
        idx = z(((k-1)*ll+1):(k*ll));
        if sum(idx==rare_state) > 0 
            w(k) = 0;
        else 
            w(k) = 1;
        end
    end
    
    for j = 1:n_latent 
        for k = 1:n_latent 
            clear wt
            if j == rare_state || k == rare_state
                wt = w;
            else 
                wt = ones(1,K);
            end
            weight_trans{j,k} = wt/sum(wt);
        end
    end
    for j = 1:n_latent 
        clear wt
        if j == rare_state
            wt = w;
        else 
            wt = ones(1,K);
        end
        weight_emit{j} = wt/sum(wt);
    end

    weights = weight(weight_trans, weight_emit);
end