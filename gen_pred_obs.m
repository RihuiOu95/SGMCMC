function y_pred = gen_pred_obs(y, trans_chain, emit_chain,  t_start, emit_dist, n_burnin, n_step_ahead, buffer, nu)
    n_mcmc = size(trans_chain.A,1);
    n_latent = size(trans_chain.A,2);
    B = buffer;
    y_pred = zeros([(n_mcmc-n_burnin),n_step_ahead]);
    z_pred = zeros([(n_mcmc-n_burnin),n_step_ahead]);
    for s = (n_burnin+1):n_mcmc
            %%%%Step 1, use forward filter to estimate p(z_T|y_1:T)
            A = reshape(trans_chain.A(s,:,:), n_latent,n_latent);
            mu = emit_chain.mu(s,:);
            sigmasq = emit_chain.sigmasq(s,:);
            if (emit_dist=="gaussian")
                emit_param = gaussian_emission_parameter(emit_chain.mu(s,:), emit_chain.sigmasq(s,:), emit_dist);
            elseif (emit_dist=="student_t")
                emit_param = student_t_emission_parameter(emit_chain.mu(s,:), emit_chain.sigmasq(s,:), nu, emit_dist);
            end
            ppred = (ones(1,n_latent)/n_latent)';
            for t = (t_start-B+1):(t_start)
                ppred = emit_param.calP(y(t))*A*ppred;
            end
            likelihood = ppred/sum(ppred);
            %%%%Step 2, simulate
            z0 = randsample(n_latent, 1, true, likelihood);
            z_pred(s-n_burnin,1) = randsample(n_latent, 1, true, A(:,z0));
            y_pred(s-n_burnin,1) = normrnd(mu(z_pred(s-n_burnin,1)), sqrt(sigmasq(z_pred(s-n_burnin,1))));
            for i = 2:n_step_ahead
                z_pred(s-n_burnin,i) = randsample(n_latent, 1, true, A(:,z_pred(s-n_burnin,i-1)));
                y_pred(s-n_burnin,i) = normrnd(mu(z_pred(s-n_burnin,i)), sqrt(sigmasq(z_pred(s-n_burnin,i))));
            end      
    end
end