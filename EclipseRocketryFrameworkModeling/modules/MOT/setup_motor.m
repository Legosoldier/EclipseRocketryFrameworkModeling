function [] = setup_motor()
%% Notes

% This script preps the motor for the Simulink block.

% 7/1/2016 -- Created

% 7/9/2016 -- Converted into function, prepping file for future creation of
% control files, building is now done solely on Simulink


%% Begin Motor Data Import

% Get list of modules used by the bessie model
%MOD = getMODstr();
%motor_name = MOD.MOT.name

% Motor selection
motor_name = 'L1520T';
motor_file = [motor_name, '.eng'];

% Import motor information
read_motor_data(motor_file);

end