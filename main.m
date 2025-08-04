% Main script for ISAC simulation
% This script reads configuration parameters and runs the ISAC simulation

clear; clc; close all;

%% Configuration
params = 'config.csv'; % Default configuration file name

% Check if configuration file exists
if ~isfile(params)
    error('Configuration file %s not found. Please create the configuration file.', config_file);
end

%% Read parameters
if ischar(params) || isstring(params)
    configTable = readtable(params);
    result = cell(height(configTable), 1);
    for i = 1:height(configTable)
        fprintf("Running configuration %d/%d...\n", i, height(configTable));
        row = configTable(i, :);
        paramStruct = table2struct(row);
        run_isac_simulation(paramStruct);
    end
    return;
end