function [transfer_rates,system_state,time] = fracflight_diffusion(end_time,input_counts,Dalpha,hfrac,alphafrac,center,grid_dim,kernex,fracorder)

% framework for a tau-leaping simulation
% from Bayati, JCP, 2013
% SSA adapted from Burrage and Carletti

global restart tauleap
global syst_type alpha_0 beta_0 rates num_channels num_nodes
global cycle_period cycle_x cycle_phase

rng('shuffle','twister');

if restart
    load temptime.mat
    current_time = time_vector(end)
else
    current_time = 0;
end

muc = 1/(2*abs(cos(alphafrac*pi/2)));

tau = 1/(12*sqrt(2));
timestep = 1;
tstepout = 1;

if strcmp(grid_dim,'1d')
    
    ffmu = (Dalpha/hfrac^alphafrac)*muc;
    
    all_a0 = zeros(1,num_nodes);
    all_propensity = zeros(num_channels,num_nodes);
    
    % fullstep refers to the moment after the diffusion takes place
    fullstep_counts = input_counts;
    % halfstep refers to the moment after tau-leaping reactions happen
    halfstep_counts = zeros(size(input_counts));
    
    transfer_rates(:,1,:) = zeros(size(input_counts));
    
    system_state(:,1,:) = input_counts;
    
elseif strcmp(grid_dim,'2d')
    
    num_nodesx = num_nodes;
    num_nodesy = num_nodes/10;
    
    input_counts = zeros(2,num_nodesx,num_nodesy)+0; %repmat(input_counts,[1,1,1,size(input_counts,3)]);
    input_counts(:,ceil(num_nodesx/2),ceil(num_nodesy/2)) = 1000;
    
    all_a0 = zeros(1,num_nodesx,num_nodesy);
    all_propensity = zeros(num_channels,num_nodesx,num_nodesy);
    
    halfstep_counts = zeros(size(input_counts));
    fullstep_counts = input_counts;
    
    hfracgrid2da = 1;
    ffmu2da = (Dalpha/hfracgrid2da^alphafrac)*muc;
    
    hfracgrid2db = sqrt(2);
    ffmu2db = (Dalpha/hfracgrid2db^alphafrac)*muc;
    
    transfer_rates(:,1,:,:) = zeros(size(input_counts));
    
    system_state(:,1,:,:) = input_counts;
    
end

while current_time < end_time
    if cycle_period>0
        alpha = alpha_0.*((1 - cycle_x) + cycle_x*sin(2*pi*current_time/cycle_period + cycle_phase));
        beta = beta_0;
        % loop over all nodes
        for idnode = 1:num_nodes
            counts = fullstep_counts(:,idnode);
            % recompute propensity vector
            all_propensity(:,idnode) = propensity(counts,alpha,beta,syst_type);
            % define basis for SSA interval
            all_a0(1,idnode) = sum(all_propensity(:,idnode));
            %         if all_a0(idnode)==0
            %             disp('there is no probability of any reaction')
            %             idnode
            %             break
            %         end
        end
    elseif strcmp(syst_type,'exponential')
        alpha = alpha_0; beta = beta_0;
        for idnode = 1:num_nodes
            counts = fullstep_counts(:,idnode);
            % recompute propensity vector
            all_propensity(:,idnode) = propensity(counts,alpha,beta,syst_type);
            % define basis for SSA interval
            all_a0(1,idnode) = sum(all_propensity(:,idnode));
            %         if all_a0(idnode)==0
            %             disp('there is no probability of any reaction')
            %             idnode
            %             break
            %         end
        end
    elseif strcmp(syst_type,'SIS')
        alpha = alpha_0; beta = beta_0;
        if strcmp(grid_dim,'1d')
            for idnode = 1:num_nodes
                counts = fullstep_counts(:,idnode);
                % recompute propensity vector
                all_propensity(:,idnode) = propensity(counts,alpha,beta,syst_type);
                % define basis for SSA interval
                all_a0(1,idnode) = sum(all_propensity(:,idnode));
                %         if all_a0(idnode)==0
                %             disp('there is no probability of any reaction')
                %             idnode
                %             break
                %         end
            end
        elseif strcmp(grid_dim,'2d')
            for idnodex = 1:num_nodes
                for idnodey = 1:num_nodes
                    counts = fullstep_counts(:,idnodex,idnodey);
                    % recompute propensity vector
                    all_propensity(:,idnodex,idnodey) = propensity(counts,alpha,beta,syst_type);
                    % define basis for SSA interval
                    all_a0(1,idnodex,idnodey) = sum(all_propensity(:,idnodex,idnodey));
                    %         if all_a0(idnode)==0
                    %             disp('there is no probability of any reaction')
                    %             idnode
                    %             break
                    %         end
                end
            end
        end
    elseif strcmp(syst_type,'noreax')
        if strcmp(grid_dim,'1d')
            all_a0 = zeros(1,num_nodes);
            all_propensity = zeros(num_channels,num_nodes);
        elseif strcmp(grid_dim,'2d')
            all_a0 = zeros(1,num_nodes,num_nodes);
            all_propensity = zeros(num_channels,num_nodes,num_nodes);
        end
    else
        if strcmp(grid_dim,'1d')
            alpha = alpha_0; beta = beta_0;
            for idnode = 1:num_nodes
                counts = fullstep_counts(:,idnode);
                % recompute propensity vector
                all_propensity(:,idnode) = propensity(counts,alpha,beta,syst_type);
                % define basis for SSA interval
                all_a0(idnode) = sum(all_propensity(:,idnode));
                %         if all_a0(idnode)==0
                %             disp('there is no probability of any reaction')
                %             idnode
                %             break
                %         end
            end
        end
    end
    
    % the SSA timestep - without the random value
    % here we need to implement the selection procedure from Cao et al,
    % which will be problematic without a further condition to consider
    % the timestep from every node, and take the minimum
    % tau leaping
    % this is a temporary and strict solution for finding a timestep that is
    % acceptable for all cities' reaction systems
    tau = min(tau,squeeze(1/max(max(all_a0))));
    tau = min(1,tau);
    %    tau = max(tau,0.01);
    
    if strcmp(grid_dim,'1d')
        
        for idnode = 1:num_nodes
            
            counts = fullstep_counts(:,idnode);
            ktau = zeros(num_channels,1);
            % scale the number of reactions during the timestep by a Poisson
            % distributed random variate
            for idx = 1:num_channels
                % Line 6 in Bayati Sec:IV.B
                poismean = tau*all_propensity(idx,idnode);
                %ktau(idx,1) = random('poiss',poismean);
                if tauleap
                    ktau(idx,1) = poissrnd(poismean);
                else
                    ktau(idx,1) = poismean;
                end
                % Line 8 in Bayati Sec:IV.B
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
        for idx = 1:num_channels
            [BigG2,tau] = frac_kernel(tau,ffmu,hfrac,num_nodes,alphafrac,center,Dalpha,fracorder);
            for lambda=1:num_nodes
                BigG2rc = BigG2(num_nodes-lambda+1:end-lambda+1);
                BigG2rc = BigG2rc./sum(BigG2rc); % renormalization: fair? it has to be normalized to work in mnrnd
                if kernex==1
                    kvector = halfstep_counts(idx,lambda).*BigG2rc; kvector = kvector'; % expectation of kernel
                else
                    kvector = mnrnd(halfstep_counts(idx,lambda),BigG2rc);
                end
                % ??? kvector(lambda) = kvector(lambda) - halfstep_counts(idx,lambda);
                vvector(idx,lambda) = vvector(idx,lambda) - sum(kvector(1:end~=lambda));
                vvector(idx,1:end~=lambda) = vvector(idx,1:end~=lambda) + kvector(1:end~=lambda);
            end
            vvector(idx,1)=0; % freeze the boundaries w.r.t. transport
            vvector(idx,end)=0;
        end
        fullstep_counts = halfstep_counts + vvector;
        %     % attempt to renormalize
        %     fullstep_counts
        %     total_pop_in
        %     sum(fullstep_counts,2)
        %     fullstep_counts = 10.*fullstep_counts;
        % step forward in time
        current_time = current_time + tau;
              
        % the index for the system trajectory
        timestep = timestep + 1;
        % keep track of the time, which is only used for diagnostics in the
        % end, since the value of end_time in ssa_flights will determine the
        % vector of time points
        time(timestep) = current_time;
        % record the system state
        %     if current_time>year+1
        %         year = year+1
        %         infecthub = tophubs(1);
        %         fullstep_counts = input_counts;
        %         halfstep_counts = zeros(num_channels,num_nodes);
        %         fullstep_counts(1,infecthub) = 0.5*input_counts(1,1,infecthub);
        %         fullstep_counts(1,infecthub) = 0.5*input_counts(1,1,infecthub);
        %     end

        if mod(round(current_time,3),1E-2)==0
            current_time
            tau
        end

        system_state(:,timestep,:) = fullstep_counts;
        transfer_rates(:,timestep,:) = vvector;
        
        %if transport=='airports'
        
        % % if transfer is largely negative, it will drive the population
        % % negative, but this must be prevented
        %     for idnode=1:num_nodes
        %         if fullstep_counts(1,idnode)<minpop;
        %             fullstep_counts(1,idnode)=minpop;
        %         end
        %         % however, the number of infected should be allowed to go to zero
        %         % if the epidemic dies away
        %         if fullstep_counts(2,idnode)<minpop;
        %             fullstep_counts(2,idnode)=0;
        %         end
        %     end
        %     min(fullstep_counts(1,:));
        
    elseif strcmp(grid_dim,'2d')
        
        % placeholder skipping reactions
        %halfstep_counts = fullstep_counts;
        
        for idx = 1:num_channels
            
            % reactions, such as SIS
            for idnodex = 1:num_nodesx
                for idnodey = 1:num_nodesy
                    
                    % this only computes infection for hubs, defined in find_probs
                    counts = fullstep_counts(:,idnodex,idnodey);
                    ktau = zeros(num_channels,1);
                    % scale the number of reactions during the timestep by a Poisson
                    % distributed random variate
                    % Line 6 in Bayati Sec:IV.B
                    poismean = tau*all_propensity(idx,idnodex,idnodey);
                    %ktau(idx,1) = random('poiss',poismean);
                    if tauleap
                        ktau(idx,1) = poissrnd(poismean);
                    else
                        ktau(idx,1) = poismean;
                    end
                    % Line 8 in Bayati Sec:IV.B
                    counts = counts + ktau(idx).*rates(:,idx);
                    if counts(idx) < 0
                        % crude way of preventing negativity - Cao et al describe a
                        % better way
                        counts(idx) = 0;
                    end
                    halfstep_counts(:,idnodex,idnodey) = counts;
                end
            end
            
            % for each node point (x0,y0), compute dispersion along four lines
            % still need to include the diagonal
            
            % fractional diffusion
            for nodepointx=1:num_nodesx % move in x, compute the y line
                [BigG2,tau] = frac_kernel(tau,ffmu2da,hfracgrid2da,num_nodesy,alphafrac,center,Dalpha,fracorder);
                XXX = nodepointx;
                lineput = squeeze(halfstep_counts(idx,XXX,:));
                vvector = zeros(num_channels,num_nodesy);
                for lambda=1:num_nodesy
                    BigG2rc = BigG2(num_nodesy-lambda+1:end-lambda+1);
                    BigG2rc = BigG2rc./sum(BigG2rc); % renormalization: fair? it has to be normalized to work in mnrnd
                    %                    kvector = mnrnd(lineput(lambda),BigG2rc);
                    kvector = lineput(lambda).*BigG2rc; kvector = kvector'; % expectation of kernel
                    % add boundary condition: rejection if hopping off grid --> if jumpsout: keep in same node end
                    % mistake in Bayati paper?kvector(lambda) = kvector(lambda) - lineput(lambda);
                    vvector(idx,lambda) = vvector(idx,lambda) - sum(kvector(1:end~=lambda));
                    vvector(idx,1:end~=lambda) = vvector(idx,1:end~=lambda) + kvector(1:end~=lambda);
                end
                fullstep_counts(idx,XXX,:) = halfstep_counts(idx,XXX,:) + reshape(vvector(idx,:),size(halfstep_counts(idx,XXX,:)));
            end
            for nodepointy=1:num_nodesy %  move in Y, compute the X line
                [BigG2, tau] = frac_kernel(tau,ffmu2da,hfracgrid2da,num_nodesx,alphafrac,center,Dalpha,fracorder);
                YYY = nodepointy;
                lineput = squeeze(fullstep_counts(idx,:,YYY)); % alternating directions
                vvector = zeros(num_channels,num_nodesx);
                for lambda=1:num_nodesx
                    BigG2rc = BigG2(num_nodesx-lambda+1:end-lambda+1);
                    BigG2rc = BigG2rc./sum(BigG2rc); % renormalization: fair? it has to be normalized to work in mnrnd
                    %                    kvector = mnrnd(lineput(lambda),BigG2rc);
                    kvector = lineput(lambda).*BigG2rc; kvector = kvector'; % expectation of kernel
                    % add boundary condition: rejection if hopping off grid --> if jumpsout: keep in same node end
                    % mistake in Bayati paper?kvector(lambda) = kvector(lambda) - lineput(lambda);
                    vvector(idx,lambda) = vvector(idx,lambda) - sum(kvector(1:end~=lambda));
                    vvector(idx,1:end~=lambda) = vvector(idx,1:end~=lambda) + kvector(1:end~=lambda);
                end
                fullstep_counts(idx,:,YYY) = fullstep_counts(idx,:,YYY) + reshape(vvector(idx,:),size(halfstep_counts(idx,:,YYY)));
            end
            
            %                         else
            %                             longpath = 1;
            %                             BigG2 = frac_kernel(tau,ffmu2da,hfracgrid2da,num_nodes,alphafrac,center,Dalpha);
            %                             if path>3
            %                                 distlB = 0;
            %                                 distrB = 0;
            %                             else
            %                                 distlB = 0;
            %                                 distrB = 0;
            %                             end
            %                         end
            
            % then compute the fractional kernel for the particular
            % line of dispersion, which is a slice through
            % halfstep_counts that is defined in a new line
            
            if path<2
            else
            end
        end
        
        %     % attempt to renormalize
        %     fullstep_counts
        %     total_pop_in
        %     sum(fullstep_counts,2)
        %     fullstep_counts = 10.*fullstep_counts;
        % step forward in time
        current_time = current_time + tau;
        
%         if mod(round(current_time,-2),0.01)==0
%             current_time
%             tau
%             time(tstepout) = current_time;
%             system_state(:,tstepout,:,:) = fullstep_counts;
%             tstepout = tstepout+1;
%         end
        
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
        %transfer_rates(:,timestep,:,:) = vvector;
        
        
    end
    
end
%size(t)
%size(X)
