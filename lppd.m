function elppd = lppd(y_rep, y, trans_chain, emit_chain,  t_start, emit_dist, breakpts, buffer, nu)
    %breakpts: the breakpoints to output log_pred_density
    n_mcmc = size(trans_chain.A,1);
    n_latent = size(trans_chain.A,2);
    n_rep = size(y_rep,1);
    n_test = size(y_rep,2);
    B = buffer;
    lpd = zeros(n_rep, n_mcmc);
    for i = 1:n_rep
        for s = 1:n_mcmc
            %%%%Step 1, use forward filter to estimate p(z_T|y_1:T)
            A = reshape(trans_chain.A(s,:,:), n_latent,n_latent);
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
            %%%%Step 2, estimate p(y_T:(T+n_test)|y_1:T)
            for t = 1 : n_test
                likelihood = emit_param.calP(y_rep(i,t))*A*likelihood;
            end
            likelihood = sum(likelihood);
            lpd(i,s) = likelihood;
        end        
    end
    
    elppd = zeros(1,length(breakpts));
    ii = 1;
    for i = breakpts
        elppd(ii) = mean(log(mean(lpd(:,1:i),2)));
        ii = ii+1;
    end
end

