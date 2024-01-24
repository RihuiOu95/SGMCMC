function [z_rep,y_rep] = gen_rep(z, n_test, n_rep,  A_true, emit_param_true)
    %%%%To generate replicate data for testing performance.
    %%%%n_rep: the number of replicate data
    z_rep = zeros([n_rep,n_test]);
    y_rep = zeros([n_rep,n_test]);
    n_latent = length(A_true);
    pi0_pred = zeros([1,n_latent]);
    pi0_pred(z(end)) = 1;
    trans_param_pred = transition_parameter(pi0_pred, A_true, A_true);
    if (emit_param_true.emit_dist == "gaussian" )
        for i = 1: n_rep
            [z_rep(i,:), y_rep(i,:)] = simHMM(trans_param_pred,emit_param_true,n_test);
        end
    elseif (emit_param_true.emit_dist == "student_t")   
        for i = 1: n_rep
            [z_rep(i,:), y_rep(i,:)] = simHMM_t_dist(trans_param_pred.pi_0, trans_param_pred.A, emit_param_true.mu, emit_param_true.sigmasq, ...
                emit_param_true.nu, n_test);
        end
    end
end