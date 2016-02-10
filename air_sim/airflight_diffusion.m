function [transfer_rates,system_state,time] = airflight_diffusion(end_time, outfreq, input_counts, transfer_prob, ...
    transfer_int, minpop, startfly_time, googleflu, flufuncarray, ffcarrayid)

% framework for a tau-leaping simulation
% from Bayati, JCP, 2013
% SSA adapted from Burrage and Carletti

global syst_type alpha_0 beta_0 rates num_channels num_nodes 
global restart tauleap
global cycle_period cycle_x cycle_phase

addpath ../viz_util/

rng('shuffle','twister');

if restart
    load temptime.mat
    current_time = time_vector(end)
else
    current_time = 0;
end

tstepout = 1;

time = zeros(size(end_time,1));

system_state(:,1,:) = input_counts;
transfer_rates(:,1,:) = zeros(size(input_counts));

% fullstep refers to the moment after the diffusion takes place
fullstep_counts = input_counts;
% halfstep refers to the moment after tau-leaping reactions happen
halfstep_counts = zeros(num_channels,num_nodes);

timestep = 1;

while current_time < end_time
    if cycle_period>0
        alpha = alpha_0.*((1 - cycle_x) + cycle_x*sin(2*pi*current_time/cycle_period + cycle_phase));
        beta_1 = beta_0;
        % loop over all nodes
        for idnode = 1:num_nodes
            counts = fullstep_counts(:,idnode);
            % recompute propensity vector
            all_propensity(:,idnode) = propensity(counts,alpha,beta_1,syst_type);
            % define basis for SSA interval
            all_a0(idnode) = sum(all_propensity(:,idnode));
            %         if all_a0(idnode)==0
            %             disp('there is no probability of any reaction')
            %             idnode
            %             break
            %         end
        end
    elseif googleflu == 1
        % loop over all nodes
        for idnode = 1:num_nodes
            if current_time > 0
                node = ffcarrayid(idnode);
                timenow = ceil(current_time*52);
                ffnow = flufuncarray(timenow,node);
                maxffnow = max(flufuncarray(:,node));
                alpha = alpha_0.*ffnow./maxffnow; % normalize rate
                beta_1 = beta_0;
                %            beta_1  = 2.5*beta_0.*((1 - .3) + 0.3*sin(2*pi*[1:1:538]/52 + 104));
                %            beta_1 = beta_0.*(1 + 2.*sin(2*pi*current_time/52/2+78).^12);
            else
                alpha = 0;
                beta_1 = 0;
            end
            counts = fullstep_counts(:,idnode);
            % recompute propensity vector
            all_propensity(:,idnode) = propensity(counts,alpha,beta_1,syst_type);
            % define basis for SSA interval
            all_a0(idnode) = sum(all_propensity(:,idnode));
            %         if all_a0(idnode)==0
            %             disp('there is no probability of any reaction')
            %             idnode
            %             break
            %         end
        end
    elseif strcmp(syst_type,'noreax')
        all_a0 = zeros(1,num_nodes);
        all_propensity = zeros(num_channels,num_nodes);
    else
        alpha = alpha_0; beta_1 = beta_0;
        for idnode = 1:num_nodes
            counts = fullstep_counts(:,idnode);
            % recompute propensity vector
            all_propensity(:,idnode) = propensity(counts,alpha,beta_1,syst_type);
            % define basis for SSA interval
            all_a0(idnode) = sum(all_propensity(:,idnode));
            %         if all_a0(idnode)==0
            %             disp('there is no probability of any reaction')
            %             idnode
            %             break
            %         end
        end
    end
    % the SSA timestep - without the random value
    % here we need to implement the selection procedure from Cao et al,
    % which will be problematic without a further condition to consider
    % the timestep from every node, and take the minimum
    % tau leaping
    % this is a temporary and strict solution for finding a timestep that is
    % acceptable for all cities' reaction systems
    tau = 1/max(all_a0);
    tau = min(transfer_int,tau);
    %tau = max(tau,0.01);
    
    for idnode = 1:num_nodes
        % this only computes infection for hubs, defined in find_probs
        counts = fullstep_counts(:,idnode);
        ktau = zeros(num_channels,1);
        % scale the number of reactions during the timestep by a Poisson
        % distributed random variate
        for idx = 1:num_channels
            % Line 6 in Bayati Sec:IV.B
            poismean = tau*all_propensity(idx,idnode);
            if tauleap
                ktau(idx,1) = poissrnd(poismean);
            else
                ktau(idx,1) = poismean;
            end            % Line 8 in Bayati Sec:IV.B
            counts = counts + ktau(idx).*rates(:,idx);
            if counts(idx) < 0
                % crude way of preventing negativity - Cao et al describe a
                % better way
                counts(idx) = 0;
            end
        end
        halfstep_counts(:,idnode) = counts;
    end
    
    vvector = zeros(num_channels,num_nodes);
    % transfer between cells/cities
    if current_time >= max(startfly_time,tau)
        if mod(current_time,transfer_int) < tau
            for idx = 1:num_channels
                for idnode = 1:num_nodes
                    kvector = squeeze(transfer_prob(idx,idnode,:)).*repmat(halfstep_counts(idx,idnode),num_nodes,1);
                    vvector(idx,idnode) = vvector(idx,idnode) - sum(kvector(1:end~=idnode));
                    vvector(idx,1:end~=idnode) = vvector(idx,1:end~=idnode) + squeeze(kvector(1:end~=idnode))';
                end
            end
        end
    end
    
    % with the ceil() the minimum transfer is unity
    fullstep_counts = halfstep_counts + sign(vvector).*abs(vvector);%ceil(abs(vvector));
    
    % if transfer is largely negative, it will drive the population
    % negative, but this must be prevented
    for idnode=1:num_nodes
        if fullstep_counts(1,idnode)<minpop;
            warning('population depleted in node %d at time %d',idnode, current_time);
            fullstep_counts(1,idnode)=minpop;
        end
        % however, the number of infected should be allowed to go to zero
        % if the epidemic dies away
        if fullstep_counts(2,idnode)<minpop;
            warning('infection depleted in node %d at time %d',idnode, current_time);
            fullstep_counts(2,idnode)=minpop;
        end
    end
    
    % step forward in time
    current_time = current_time + tau;
    
    if mod(timestep,outfreq)==0
        time(tstepout) = current_time;
        system_state(:,tstepout,:) = fullstep_counts;
        transfer_rates(:,tstepout,:) = vvector;
        tstepout = tstepout+1;
    end
    if mod(timestep,20)==0
        current_time
    end
    % the index for the system trajectory
    timestep = timestep + 1;
    % keep track of the time, which is only used for diagnostics in the
    % end, since the value of end_time in ssa_flights will determine the
    % vector of time points
    % record the system state
    %     if current_time>year+1
    %         year = year+1
    %         infecthub = tophubs(1);
    %         fullstep_counts = input_counts;
    %         halfstep_counts = zeros(num_channels,num_nodes);
    %         fullstep_counts(1,infecthub) = 0.5*input_counts(1,1,infecthub);
    %         fullstep_counts(1,infecthub) = 0.5*input_counts(1,1,infecthub);
    %     end
end
%size(t)
%size(X)
