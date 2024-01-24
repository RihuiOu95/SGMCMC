z = [];
for i=1:length(data_cell)
    if i==4
       z = [z; impute(data_cell{1,i},data_cell{1,i})];
    else
       z = [z; impute(data_cell{1,i},data_cell{2,i})];
    end
end
%%% 5-moving average
z = z(19:end);
smoothz = mean(reshape(z,[1,length(z)/1]),1);
logz = log10(smoothz);
figure
plot(logz)
title('Original Data')
