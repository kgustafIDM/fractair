clear all

nodeid= 0; % node to be infected for delta-function outbreak
platform = 2; % unix/mac == 1, windows == 2

if platform==1
    dumfile = '~/Dropbox/EMOD/results/filefile.fl';
elseif platform==2
    dumfile = '.\filefile.fl';
end
num_nodes = 482;
for nodeid = 1:num_nodes;
    runname = air_paramscan23(nodeid,platform)
    while exist(dumfile,'file')
        pause(1)
    end
end