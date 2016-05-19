function [googleflu,country,startfly_time,outfreq,transfer_int,num_nodes,total_indiv,transfer_prob,gflu51states,idState,nodeidUSA,nodecodesUSA] = airflight_input

% input parameters for air network diffusion
% meant to be modified by user for each production run

% optional: use the open source fftf filtering program to get top 10 Fourier coeffs
% by Fourier transforming the Google Flu data for the entire country

% to activate the Google flu data
googleflu = 1; % on is 1

% many countries are available, please inquire to kgustafson@idmod.org to discuss
country = 'USA'
% placeholder
infecthub = 0;
gflu51states = 1; 
idState = 0;
popbase = 1000;

% time at which to start the transport
startfly_time = 0;
% frequency of outputs
outfreq = 1;
% interval between transfers due to air travel
transfer_int = 0.0192; %0.0192 ???
% only a single city ?
singlecity = 0;
% for a fake simulation ?
fakey=0;

% probability matrix of transfer from city to city (cell to cell)
% how to determine the probability of remaining?
% v_lambda; v_lambda-1; v_lambda+1
% these test cases are set up so that the row indicates output travelers
% from that city to the nth city

% this factor is scaling the transfer matrix to account for the rescaling
% of the population of the airport service regions - but the interval
% between transports is equally important
% see tauleapdiffusion.m
trfactor = 1.5/1E7; 

%trfactor 1/1E7; works for USA noreax at 7 nodes for a short time

% probability matrix of transfer from city to city (cell to cell)
% how to determine the probability of remaining?
% v_lambda; v_lambda-1; v_lambda+1
% these test cases are set up so that the row indicates output travelers
% from that city to the nth city


if fakey==1
    num_nodes = 4;
    transfer_prob1 = [0 0; 2 2; 2 2; 2 2];
    transfer_prob2 = [0 0; 0 0; 0 0; 0 0];
    transfer_prob3 = [0 0; 0 0; 0 0; 0 0];
    transfer_prob4 = [4 4; 0 0; 0 0; 0 0];
    
    transfer_prob = cat(3,transfer_prob1',transfer_prob2',transfer_prob3',transfer_prob4').*trfactor;
else
    if strcmp(country,'USA') == 1
        % using the 482 nodes with more than one connection
        load ../viz_util/flights_network/countries/USA/data/transferUSAnodes.mat
        % using the top 7 nodes
        %load /Users/gustafsonkyle/Research/EMOD/Basil1/flights_network/countries/USA/data/transferUSAnodes7.mat
        load ../viz_util/flights_network/countries/USA/data/gflu51states.mat
        % symmetrical transfer just corrects for small inconsistencies in
        % the OAG flight data
        transfer_prob = symtransferUSA.*trfactor;
        % infected airport city ids
        infecthub = tophubsUSA(1);
        % count_unique is an add-on function
        num_nodes = size(nodeidUSA,1);
        nn = origincount(nodeidUSA);
        total_indiv = nn + popbase; %round(1000.*nn./mean(nn(tophubsUSA)))+100; % 100 sets a minimum population size
    end
end

if singlecity==1
    transfer_prob = transfer_prob(:,:,infecthub);
end
