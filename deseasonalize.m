
function [res] = deseasonalize(y, period)
    T = length(y);
    sW5000 = repmat(1/period,period,1);
    yS = conv(y,sW5000,'same');
    yS(1:period) = yS(period+1); yS(T-period:T) = yS(T-period-1);
    xt = y-yS;
    %%%%%
    s = period-1;
    sidx = cell(s,1);
    for i = 1:s
     sidx{i,1} = i:s:T;
    end
    sst = cellfun(@(x) mean(xt(x)),sidx);
    % Put smoothed values back into a vector of length N
    nc = floor(T/s); % no. complete years
    rm = mod(T,s); % no. extra months
    sst = [repmat(sst,nc,1);sst(1:rm)];

    % Center the seasonal estimate (additive)
    sBar = mean(sst); % for centering
    sst = sst-sBar;
    %%%%
    res = y - sst;
end