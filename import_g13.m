%% Import data from text file
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