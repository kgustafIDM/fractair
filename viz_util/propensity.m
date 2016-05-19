function a = propensity(input_counts,alpha,beta,syst_type)

% returns propensities for present state
% species, from Schmelzer, Burrage, Carletti

global species1 species2 

total_indiv = sum(input_counts);
    
if strcmp(syst_type,'SIR')
           
    if total_indiv > 0 
        c = [alpha/total_indiv beta];
    else
        c = [0 beta];
    end
    
    input_counts = input_counts';
    
    a(1) = c(1)*input_counts(species1)*input_counts(species2);
    a(2) = c(2)*input_counts(species2);
    
    a=[a(1) a(2)]';
elseif strcmp(syst_type,'SI')
           
    if total_indiv > 0 
        c = [alpha/total_indiv beta];
    else
        c = [0 beta];
    end
    
    input_counts = input_counts';
    
    a(1) = c(1)*input_counts(species1)*input_counts(species2);
    a(2) = c(2)*input_counts(species2);
    
    a=[a(1) a(2)]';
    
elseif strcmp(syst_type,'SIS')
    % apparently the same propensity vector can be used for SIS as for SIR
    % just need to change the stochiometry in SIR_flights.m
                   
    if total_indiv > 0 
        c = [alpha/total_indiv beta];
    else
        c = [1 beta];
    end
    
    input_counts = input_counts';
    
    a(1) = c(1)*input_counts(species1)*input_counts(species2);
    a(2) = c(2)*input_counts(species2);
    
    a=[a(1) a(2)]';
elseif strcmp(syst_type,'exponential')
                
    if total_indiv > 0 
        c = [alpha/total_indiv]; % where alpha is rho for growth rate
    else
        c = [0];
    end
    
    input_counts = input_counts';
    
    a(1) = c(1)*input_counts(species1);
    a(2) = c(1)*input_counts(species2);
    
    a=[a(1) a(2)]';
else
    error('invalid system type')
end