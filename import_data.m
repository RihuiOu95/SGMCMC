
%% Import g13 data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/User/MATLAB/projects/sgmcmc_new/OOP_matlab/codes/201709_g13.csv
%
% Auto-generated by MATLAB on 03-Sep-2020 22:41:25

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [169, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["time_tag", "A_QUAL_FLAG", "A_NUM_PTS", "A_AVG", "B_QUAL_FLAG", "B_NUM_PTS", "B_AVG"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
data09_g13 = readtable("data/201709_g13.csv", opts);
data10_g13 = readtable("data/201710_g13.csv", opts);
data11_g13 = readtable("data/201711_g13.csv", opts);
data12_g13 = readtable("data/201712_g13.csv", opts);
data01_g13 = readtable("data/201801_g13.csv", opts);
data02_g13 = readtable("data/201802_g13.csv", opts);
data03_g13 = readtable("data/201803_g13.csv", opts);
data04_g13 = readtable("data/201804_g13.csv", opts);
%% Clear temporary variables
clear opts



%% Import g15 data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/User/MATLAB/projects/sgmcmc_new/OOP_matlab/codes/201709.csv
%
% Auto-generated by MATLAB on 27-Aug-2020 19:09:23

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [168, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["time_tag", "A_QUAL_FLAG", "A_NUM_PTS", "A_AVG", "B_QUAL_FLAG", "B_NUM_PTS", "B_AVG"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
data09 = readtable("data/201709.csv", opts);
data10 = readtable("data/201710.csv", opts);
data11 = readtable("data/201711.csv", opts);
data12 = readtable("data/201712.csv", opts);
opts.DataLines = [169, Inf];
data01 = readtable("data/201801.csv", opts);
data02 = readtable("data/201802.csv", opts);
data03 = readtable("data/201803.csv", opts);
data04 = readtable("data/201804.csv", opts);
%%%Store in the data cell
indices = ["09","10", "11", '12', '01', '02', '03', '04'];
data_cell = cell(2, length(indices));
for i = 1 : length(indices)
    index = indices(i);
    data_cell{1,i} = eval(strcat("data",index,".B_AVG"));
    data_cell{2,i} = eval(strcat("data",index,"_g13",".B_AVG"));
end
% Delete Missing data

%y = y(136:end);
%smoothy = mean(reshape(y,[5,length(y)/5]),1);
%% Clear temporary variable
clear opts