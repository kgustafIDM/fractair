function [grid_dim,Dalpha,hfrac,center,alphafrac,num_nodes,total_indiv,infecthub,step_infect,kernex,fracorder] = fracflight_input

% input parameters for fractional diffusion

grid_dim = '1d';
hfrac = 0.01;

%grid_dim = '2d';

Dalpha = 1;
center = 1;
alphafrac = 0.3;
fracorder = 2;
num_nodes = 482;
total_indiv = repmat(1000,1,num_nodes);

% temporary num_nodes = size(nodecodesUSA,1);
infecthub = ceil(num_nodes/2);
step_infect = 0; 
kernex = 1; % yes/no using expectation of fractional kernel instead of multinomial sampling
