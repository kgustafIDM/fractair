% runairparamscan.m will run several instances of the air flight compartmental (SIS/SIR) model
% Most of the parameters for the compartmental model are contained in air_paramscan.m or downstream files
% Dependencies include the pre-compiled transfer matrix for the air flight database in the country of interest

clear all

nodeid= 0; % node to be infected for delta-function outbreak
platform = 2; % unix/mac == 1, windows == 2

if platform==1
    % this is a dummy file as a placeholder to force a pause the prevents a crash
    dumfile = '~/Dropbox/EMOD/results/filefile.fl';
elseif platform==2
    dumfile = '.\filefile.fl';
end

% the number of nodes is actually the set of nodes for testing the initial position of the outbreak infection
% there are several other parameters for controlling the position and amplitude of the outbreak node, these can be found in downstream files
infect_nodes = 1:482;
for nodeid = infect_nodes;
    runname = air_paramscan23(nodeid,platform)
    while exist(dumfile,'file')
        pause(1)
    end
end
