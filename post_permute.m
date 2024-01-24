function [A_sorted, mu_sorted, sigmasq_sorted] = post_permute(emit_chain,trans_chain)
    %permute posterior samples by mu_1<mu_2<...<mu_K
    n_mcmc = size(emit_chain.mu,1);
    n_latent = size(emit_chain.mu,2);
    sigmasq_sorted = zeros(n_mcmc,n_latent);
    A_sorted = zeros(n_mcmc, n_latent,n_latent);
    %%%
    [mu_sorted, idx] = sort(emit_chain.mu,2);
    for i = 1:n_mcmc
        sigmasq_sorted(i,:) = emit_chain.sigmasq(i,idx(i,:));
        A = reshape(trans_chain.A(i,:,:),n_latent,n_latent);
        A_sorted(i,:,:) = A(idx(i,:),idx(i,:));
    end
end