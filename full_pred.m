function [pipred,qpred] = full_pred(y, t, trans_param, emit_param)
    T = size(y,1);
    n_latent = emit_param.n_latent;
    A = trans_param.A;   
    qpred = eye(n_latent);
    pipred = trans_param.pi_0';
    if t == 1
        for i = (t+1):T
            qpred = emit_param.calP(y(i))*A*qpred;
        end
    elseif t == T
        for i = (t+1):T
            qpred = emit_param.calP(y(i))*A*qpred;
        end
    else
        for i = 1:(t-1)
            pipred = emit_param.calP(y(i))*A*pipred;
        end
        for i = (t+1):T
            qpred = emit_param.calP(y(i))*A*qpred;
        end
    end
    qpred = ones(n_latent,1)'*qpred;
    qpred = qpred';
end