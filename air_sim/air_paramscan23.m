function runname = air_paramscan(nodeid,platform)

% platform: unix/mac=1, windows=1 

if platform==1
    !touch ../results/filefile.fl
elseif platform==2
    !type NUL > ..\results\filefile.fl
end

% Several important parameters are set to be global for convenience and
% safety.

global restart tauleap
global species1 species2 syst_type alpha_0 beta_0 rates num_channels num_nodes
global init_inf_percent total_indiv
global cycle_period cycle_x cycle_phase

% The nodeid is typically specified in runairparamscan.m as part of a
% parameter scan through different initial infected nodes. Of course this
% is a parameter scan only in the sense that the initial location of the
% infection is a parameter.

infecthub = nodeid;

% Another possible initial condition for the infection, useful for studying
% the wavefront of the infection propagation, is a step-wise infection. The
% step is a bit difficult to define on the air travel network, so it has
% been hardwired below as an infection in the first 1/3 of the list of
% nodes, which is sorted alphabetically, so this is not an actual
% geographic step, but rather a curiosity. 

% this should usually be set to zero
step_infect = 0;


% It is possible to restart a paused simulation if an intermediate
% simulation state has been saved as "temptime.mat". This will cause
% intermediate input_counts and an intermediate time to be set as the
% initial conditions.

restart = 0; % if set to 1 (on), then a temptime.mat input file must be in the running folder

% Here, an important choice is made between doing a stochastic simulation,
% where a Poisson random number is chosen around the mean value for a
% reaction propensity, or if the mean value is used instead. See
% airflight_diffusion and Bayati JCP 2013 for more details.

tauleap = 0; % if set to 1 (on), the reaction process will be stochastic

% total time in years
end_time = 5;

% Choose the type of compartmental model for the disease: SI/SIR/SIS are
% supported and specified through the stoichiometric rate matrices found
% below.

syst_type = 'SIS';

% Cycling parameters for seasonality, which is implemented as an
% oscillation in alpha within the airflight_diffusion backend. 

cycle_period = 0; % units of years
cycle_x = 0.5;
cycle_phase = 0.0;

% Set the important parameters for the SIR or SIS compartmental model
% notice that the R_0 value is tunable, rather than the alpha_0 

R_0 = 2;  % reproductive number magnitude
beta_0 = 3; % recovery rate
alpha_0 = R_0*beta_0; % infection rate per unit infected

% Compartmental model rate matrix (stoichiometry)
if strcmp(syst_type,'SI')
    rates = [-1 0; 1 0];
elseif strcmp(syst_type,'SIR')
    rates = [-1 0; 1 -1];
elseif strcmp(syst_type,'SIS')
    % a symmetric, sometimes oscillatory system
    rates = [-1 1; 1 -1];
elseif strcmp(syst_type,'noreax')
    rates = [0 0; 0 0];
end

species1 = 1; % susceptible
species2 = 2; % infected

[num_species,num_channels] = size(rates);

% use the open source fftf filtering program to get top 10 Fourier coeffs
% for the google flu data for the entire country

% Initialize the infection, defaults

init_inf_percent = 50; % measure infected percentage-wise
total_indiv = 1000; %round(1000.*nn./mean(nn(tophubsUSA)))+100; % 100 sets a minimum
%infecthub = 1; % default infected hub

% the deltafect parameter should be ignored
deltafect = 1; % all other hubs have zero infection
% requiring a floor on population
minpop = 0;

%[grid_dim,dumdum,center,bumbum,num_nodes,infecthub,step_infect,kernex] = fracflight_input;

% Call the file with most of the input parameters

[googleflu,~,startfly_time,outfreq,transfer_int,num_nodes,total_indiv,transfer_prob,gflu51states,idState,nodeidUSA,nodecodesUSA] = airflight_input;

% initial state baseline
initial_states = zeros(num_species,num_nodes);
% the total number of individuals using the hub can be scaled according to the
% connectivity of the hub
init_inf = round(total_indiv*init_inf_percent/100); % set the number of infected
%initial_states(species1,:) = total_indiv-init_inf;
%initial_states(species2,:) = init_inf; % assign number of infected
% a special infected hub can be selected
if infecthub>0
    if deltafect>0
        initial_states(1,:) = total_indiv;
        initial_states(2,:) = 0;
        initial_states(1,infecthub) = total_indiv(infecthub)-init_inf(infecthub);
        initial_states(2,infecthub) = init_inf(infecthub);
    else
        initial_states(:,infecthub) = [0.5*sum(initial_states(:,infecthub),1);...
            0.5*sum(initial_states(:,infecthub),1)];
    end
elseif infecthub==0 && step_infect==1
    initial_states(1,1:floor(size(initial_states,2)/3)) = total_indiv-init_inf;
    initial_states(2,1:floor(size(initial_states,2)/3)) = init_inf;
    initial_states(1,1+floor(size(initial_states,2)/3):end) = total_indiv;
    initial_states(2,1+floor(size(initial_states,2)/3):end) = 0;
end

[transfer_rates,final_state,time_vector] = airflight_diffusion(end_time,outfreq,initial_states,transfer_prob,transfer_int,minpop,startfly_time,googleflu,gflu51states,idState);
runname = sprintf('air23_%d',nodeid);
timev = time_vector;
timev = timev(1:20:end);
finals = final_state;
finals = finals(2,1:20:end,:);
meanfinalnorm = squeeze(mean(final_state(2,:,:),3))./1000;
h1=figure('Visible','off'); plot(time_vector,meanfinalnorm,'k','Displayname',runname);
xlabel('time (yr)'); ylabel('fraction infected'); title(nodecodesUSA{nodeidUSA(nodeid)});
%sdf('ploscompbio');
hold on; legend('show');

% sigwork = squeeze(mean(final_state(2,:,:),3))./1000;
% sizesigwork = size(sigwork,2);
% intwork = trapz(time_vector,sigwork);
% clear final_state time_vector

% if platform==1
%     load ~/Research/EMOD/SISfrac/SISUSA_2_3_nocycle.mat final_state time_vector
% elseif platform==2
%     load C:\Users\kgustafson\Dropbox\EMOD\results\SISfrac\SISUSA_2_3_nocycle.mat final_state time_vector
% end

% find index of the USA data for the end_time specified in the scan
% endtimeid = find(min(abs(time_vector-end_time))==abs(time_vector-end_time));

% plot(time_vector(1:endtimeid),squeeze(mean(final_state(2,1:endtimeid,:),3))./1000,'r','DisplayName','USA (2,3)');
saveas(h1, [runname,'.png']);
saveas(h1, [runname,'.fig']);
close(h1);

% sigUSA23 = squeeze(mean(final_state(2,1:endtimeid,:),3))./1000;
% sizesigUSA23 = size(sigUSA23,2);
% intUSA23 = trapz(time_vector(1:endtimeid),sigUSA23);

%diffworkUSA = intwork-intUSA23;

% necessary to resample signals for computing the distance
% if sizesigUSA23<sizesigwork
%     resamp = floor(linspace(1,sizesigwork,sizesigUSA23));
%     sigwork = sigwork(1,resamp);
% elseif sizesigUSA23>sizesigwork
%     resamp = floor(linspace(1,sizesigUSA23,sizesigwork));
%     sigUSA23 = sigUSA23(1,resamp);
% else
%     sigUSA23 = sigUSA23;
%     sigwork = sigwork;
% end
% 
% eucdistUSA = pdist([sigUSA23;sigwork]);
% l2normUSA = cumsum((sigwork-sigUSA23).^2);
save([runname,'_gflu.mat'],'init_inf_percent','deltafect','meanfinalnorm','finals','timev','googleflu','num_nodes','syst_type','R_0','beta_0');

if platform==1
    !rm ../results/filefile.fl
elseif platform==2
    !del ..\results\filefile.fl
end