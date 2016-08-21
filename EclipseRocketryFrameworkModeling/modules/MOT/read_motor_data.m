function MOT_data = read_motor_data(fileName)
% Reads .eng files and pulls out thrust data vs. time, as well as some
% metadata about the motor and returns the info as a struct. 
% Motor data files come from rocketreviews.com

% check if fileName exists
if ~exist(fileName, 'file')
    error('File does not exist!!\n');
    return
end

% check if fileName is proper input
if fileName(end-3:end) ~= '.eng'
    error('File is not a .eng file!\n');
    return
end

% read file into MATLAB
file_dat = textread(fileName, '%s', 'delimiter', '\n');

% extract meta-data from file
metaDat = strsplit(file_dat{3});
name = metaDat{1};
diameter = str2double(metaDat{2}); % mm
length = str2double(metaDat{3}); % mm
prop_weight = str2double(metaDat{5}); % grams 
tot_weight = str2double(metaDat{6}); % grams
mfger = metaDat{7};

% extract time and thrust data
time = [];
thrust = [];
nRows = numel(textread(fileName,'%1c%*[^\n]'));
for i = 4:(nRows-1)
   line = strsplit(file_dat{i});
   time(end+1) = str2double(line{1});
   thrust(end+1) = str2double(line{2});
end

% calculate the total impulse
impulse = spline(time, thrust);
tot_imp = integral(@(x)ppval(impulse,x),time(1), time(end));

% initialize and pass data into struct
MOT_data = struct('name', name,...
                    'mfger', mfger,...
                    'time', time,...
                    'thrust', thrust,...
                    'diameter', diameter,...
                    'length',length,...
                    'prop_wt', prop_weight,...
                    'tot_imp', tot_imp,...
                    'tot_wt',tot_weight);

% save into .mat file for later useage
save MOT_lib.mat MOT_data
end