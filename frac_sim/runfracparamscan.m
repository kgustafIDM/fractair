clear all

Dalpha = 0;
alphafrac = 0;

platform = 2; % unix/mac == 1, windows == 2

if platform==1
    dumfile = '~/Dropbox/EMOD/results/filefile.fl';
    addpath ../functions/
elseif platform==2
    dumfile = '.\filefile.fl';
    addpath ..\viz_util\
end

for Dalpha = [2]
    for alphafrac = 0.1:0.1:0.2
        Dalpha
        alphafrac %#ok<*NOPTS>
        runname = frac_paramscan(Dalpha,alphafrac,platform)
        while exist(dumfile,'file')
            pause(1)
        end
    end
end