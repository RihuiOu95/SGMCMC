function logpredden = comp_logpredden(alt_data_y, pi0true, Asave, musave, sigmasqsave, ttest, S)
    %Calculate p(y_{T:(T+ttest)}|y_1:T) = \sum_{z_T} p(y_{T:(T+ttest)}|z_T)*p(z_T|y_{1:T})
    % Output: the expected log predictive density
    % Input: alt_data_y: ytest 
    % pi0_true: the K by 1 distribution vector p(z_T|y_{1:T}) which can be calculated by forward filter
    % Asave: the posterior samples of A
    % musave: the posterior samples of mu  
    % sigmasqsave: the posterior samples of sigmasq
    % ttest: the length of a test data, set to be 10 here
    % S: the number of posterior samples used in this case. We only use a subsample because of effficieny
    
    post_size = length(musave); %posterior sample size
    n = length(alt_data_y); % the numbers of alternative test data; 
    % note that it is different from the length of a test data here
    logpredden = 0; %Initialize the logpredictive density variable
        for i = 1:n
        predden=zeros(S,1);
        schoose = randsample(post_size,S); %Use a subset of posterior samples to increase efficieny
        ii = 1;
        for s = schoose'
            likelihood = pi0true'; %Initialize the predictive likelihood by p(z_T|y_{1:T})
            for t = 1:ttest
                likelihood = calP(alt_data_y(i,:),sigmasqsave(s,:),musave(s,:),t)*Asave(:,:,s)*likelihood;
                %Use the foruma of p(y_{T:{T:T+ttest}}|z_T)=\Pi_{t=T+1}^{t=T+ttest}P(y_t)A
            end
            likelihood = sum(likelihood); %Marginalizing out z_T
            predden(ii) = likelihood;
            ii = ii+1;
        end
        predden = mean(predden); %calculate the mean of predictive density, i.e., sum_{s=1}^S p(y|theta^s)/S
        logpredden = logpredden + log(predden);
    end 
end