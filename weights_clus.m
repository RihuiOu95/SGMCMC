function weights = weights_clus(y, n_latent, K, L, rare_state, eta) 

    [index, center] = kmeans(y, n_latent, 'MaxIter', 1000); %Kmeans with K clusters in total
    clear n weights_trans wt_emit weights_emit
    for i = 1:n_latent
        n(i) = sum(index==i); 
    end
    
    for k = 1:K
        idx = index(((k-1)*(2*L+1)+1):(k*(2*L+1)));
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