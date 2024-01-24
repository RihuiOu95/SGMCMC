function minibatch = form_mb(clusters, b)
    clear mb;
    for j = 1:clusters.n_clusters
        mb{j} = randsample(clusters.mb{j},b(j));
    end
    minibatch = mini_batch(mb);
end