function runname = frac_paramscan43(Dalpha,alphafrac,platform)

% platform: unix/mac=1, windows=2 

% simple reaction simulation
% from Bayati JCP 2013

% 1.) a = (a1,a2,...)^T : probability of reactions
% 2.) deltat_r = Tt : sample from exponential
% 3.) l = Qq ~ a_j(t)/a_0(t) : sample point-wise distribution of a
% 4.) n_bar(t + deltat_r) = n_bar(t) + nu_l

% multinomial stochastic algorithm with Lie-Trotter splitting and tau-leap

% 1.) a = (a1,a2,...)^T : the probablility of reactions
% 2.) deltat_r = Tt : sample from exponential or other type of reaction
% 3.) deltat_d = Dd : time-step restriction
% 4.) deltat = min(deltat_r,deltat_d)
% for all i do
% 5.) k_i ~ Pp(a_i*deltat) : number of executions for each channel i
% end for
% 6.) n_bar(t + deltat/2) = n_bar(t) + sum_l[k_l*nu_l] : state vector,
%                                                        number of particles
% 7.) nu_d = 0 : diffusion transitions
% for all lambda do : loop over cells, spatial discretization
% 8.) k1 ~ Bb(n_barlambda(t + deltat/2),Gg_lambdam1_lambda(deltat) :
%           binomial distribution Bb(L,p) of L trials, probability p
% 9.) k2 ~ Bb(n_barlambda(t + deltat/2) - k1,
%               (Gg_lambdap1_lambda(deltat))/(1-Gg_lambdam1_lambda(deltat))
% 10.) nu_lambda_d = nu_lambda_d - (k1 + k2) : movement out of lambda
% 11.) nu_lambdam1_d = nu_lambdam1_d + k1 : movement into neighbor
% 12.) nu_lambdap1_d = nu_lambdap1_d + k2 : movement into other neighbor
% end for
% 13.) n_bar(t + deltat) = n_bar(t + deltat/2) + nu_d : new state vector

if platform==1
    !touch ~/Dropbox/EMOD/results/filefile.fl
elseif platform==2
    !type NUL > .\filefile.fl
end

global restart tauleap
global species1 species2 syst_type alpha_0 beta_0 rates num_channels num_nodes
global init_inf_percent total_indiv
global cycle_period cycle_x cycle_phase

grid_dim = '1d';
%grid_dim = '2d';

center = 1;

% temporary num_nodes = size(nodecodesUSA,1);
num_nodes = 482;
infecthub = ceil(num_nodes/2);
step_infect = 0;
kernex = 1; % yes(1)/no(2) using expectation of fractional kernel instead of multinomial sampling

hfrac = 1; 
fracorder = 2;

restart = 0;

tauleap = 1;

% total time in years
end_time = 1;

% input parameters for reactions

syst_type = 'SIS';

% allow alpha to cycle for seasonality
cycle_period = 0; % units of years
cycle_x = 0.5;
cycle_phase = 0.0;

R_0 = 4;  % reproductive number magnitude
beta_0 = 3; % recovery rate
alpha_0 = R_0*beta_0; % infection rate per unit infected

% rate matrix (stoichiometry)
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

init_inf_percent = 100; % measure infected percentage-wise
total_indiv = 1000; %round(1000.*nn./mean(nn(tophubsUSA)))+100; % 100 sets a minimum
%infecthub = 1; % default infected hub
deltafect = 1; % all other hubs have zero infection
% requiring a floor on population
minpop = 0;

%[grid_dim,dumdum,center,bumbum,num_nodes,infecthub,step_infect,kernex] = fracflight_input;

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
        initial_states(1,infecthub) = total_indiv-init_inf;
        initial_states(2,infecthub) = init_inf;
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


[transfer_rates,final_state,time_vector] = fracflight_diffusion(end_time,initial_states,Dalpha,hfrac,alphafrac,center,grid_dim,kernex,fracorder);
runname = sprintf('frac_fig2noise_%.2f_%.2f',Dalpha,alphafrac);
%timev = time_vector; 
% 
% h1=figure('Visible','off'); plot(time_vector,squeeze(mean(final_state(2,:,:),3))./1000,'k','Displayname',runname);
% hold on; legend('show');
% 
% sigwork = squeeze(mean(final_state(2,:,:),3))./1000;
% sizesigwork = size(sigwork,2);
% intwork = trapz(time_vector,sigwork);
% clear final_state time_vector
% 
% if platform==1
%     load ~/Research/paperBayati/SISfrac/SISUSA_4_3_nocycle.mat final_state time_vector
% elseif platform==2
%     load C:\Users\kgustafson\Dropbox\paperBayati\results\SISfrac\SISUSA_4_3_nocycle.mat final_state time_vector
% end
% 
% % find index of the USA data for the end_time specified in the scan
% endtimeid = find(min(abs(time_vector-end_time))==abs(time_vector-end_time));
% 
% plot(time_vector(1:endtimeid),squeeze(mean(final_state(2,1:endtimeid,:),3))./1000,'r','DisplayName','USA (2,3)');
% saveas(h1, [runname,'.png']);
% close(h1);
% 
% sigUSA43 = squeeze(mean(final_state(2,1:endtimeid,:),3))./1000;
% sizesigUSA43 = size(sigUSA43,2);
% intUSA43 = trapz(time_vector(1:endtimeid),sigUSA43);
% 
% diffworkUSA = intwork-intUSA43;
% 
% % necessary to resample signals for computing the distance
% if sizesigUSA43<sizesigwork
%     resamp = floor(linspace(1,sizesigwork,sizesigUSA43));
%     sigwork = sigwork(1,resamp);
% elseif sizesigUSA43>sizesigwork
%     resamp = floor(linspace(1,sizesigUSA43,sizesigwork));
%     sigUSA43 = sigUSA43(1,resamp);
% else
%     sigUSA43 = sigUSA43;
%     sigwork = sigwork;
% end
% 
% eucdistUSA = pdist([sigUSA43;sigwork]);
% l2normUSA = cumsum((sigwork-sigUSA43).^2);
% save([runname,'.mat'],'finals','timev','sigwork','l2normUSA','eucdistUSA','diffworkUSA','intwork','intUSA23','Dalpha','alphafrac','num_nodes','syst_type','R_0','beta_0');

save([runname,'.mat'],'time_vector','final_state','Dalpha','alphafrac','num_nodes','syst_type','R_0','beta_0','kernex','tauleap');

if platform==1
    !rm ./filefile.fl
elseif platform==2
    !del .\filefile.fl
end