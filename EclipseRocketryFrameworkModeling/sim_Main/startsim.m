% setups and runs the whole simulation
% provides the parameters for the model

clear;
clc;
close all;
bdclose all;

setup_filepaths;

%% INPUT PARAMETERS
% Here, we input the scenarios to modify the simulation parameters
motor_name = 'L1520';
% scenario = motor_name;

%% Load input parameters

% get relevant modules
%MODs = getMODstr();

% setup all modules
%setup_all_mods;

% pass input parameters into workspace
%load_all_mods;

%% Run Simulation

% Open up simulink model
sim_Name = 'Bessie';
sim_File = ([sim_Name '_Main']);
open_system(sim_File);

% run simulation
sim_Out = sim(sim_File, 'SimulationMode', 'normal');

%% Plot Important Results

% figure
% plot stuffhere
