function [] = setup_filepaths()
% this function will load all the necessary files into MATLAB's search
% directory

% insert section for gathering the proper modules we will need

% Get directory data
[dir, ~,~] = fileparts(pwd);

% load important file locations aside from modules needed
addpath(genpath(dir));

% change to main directory
cd ../

end