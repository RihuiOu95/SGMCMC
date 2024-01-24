function stat = cal_stationary(trans_chain,iter_num)
    n_latent = length(trans_chain.A(iter_num,:,:));
    A = reshape(trans_chain.A(iter_num,:,:),n_latent,n_latent);
    [u,v]=eig(A);
    index = find(abs(diag(v)-1)<1e-7);
    stat = u(:,index) ./ sum(u(:,index));
end