function [len, colony_history, nutrient_history, absorb_history,...
    final_path] = nutrient_model_colony(init_size, time, ...
    c1_concentration, move_prob, absorb_prob, replicate_prob,...
    pathline, pathline_vals)

% == NUTRIENT MODEL COLONY ==
%   Model of cylindrical yeast growth with nutrient diffusion.
%
% = INPUTS =
%   o init_size: initial size of the yeast colony
%   o time: number of discrete time periods for simulation
%   o pathline: boolean indicating whether to graph pathline plots
%   o pathline_vals: values of cells to create pathline plots for
%   o c1_concentration: concentration of nutrient in first cell
%   o move_prob: chance of nutrient movement
%   o absorb_prob: chance of nutrient absorption by yeast
%   o replicate_prob: chance of proliferation per nutrient

% = OUTPUTS =
%   o len: length of the colony at each time period
%   o colony_history: history of the colony at each time period with cells
%   named
%   o nutrient_history: available nutrient in each position at each time period
%   o absorb_history: nutrient absorbed by each cell at each time period
%   o final_path: pathline of selected cells
%       
% See also: nmc_repetitions, nutrient_distribution



% ======================================================================
% Default values for pathline and pathline_vals
% ======================================================================

switch nargin
    case 3
        pathline = false;
        pathline_vals = 1;
    case 4
        pathline_vals = 1;
end

% ======================================================================
% Initialising variables
% ======================================================================

% Array of numbered yeast cells, up to given size
colony = 1:init_size;

% Storage for colony data at each time step, including at t = 0
colony_history = cell(time+1, 1);
colony_history{1} = colony;

% (Add description here)
nutrient_history = cell(time+1, 1);
absorb_history = cell(time+1, 1);

% Next value to be added to colony
next_cell = init_size + 1;

% Storage for length data at each time step, including at t = 0
len = init_size;

% Define path if requested
if pathline
    
    % Initialise as zero matrix:
    % row = position of a yeast cell
    % col = time step
    path = zeros(size(pathline_vals,2),time+1);
    for i = 1:size(pathline_vals,2)
        
        % Do not include position of values not initially in colony
        if pathline_vals(i)>init_size
            path(i,1) = NaN;
        
        % Find initial position of values requested
        else
            path(i,1) = find(colony == pathline_vals(i));
        end
    end
    
else
    final_path = 0;
end

% Initialise nutrient distribution variables
nutri_c1_concentration = c1_concentration;
nutri_movement_prob = move_prob;
nutrient = zeros(1,init_size);
nutrient(1) = nutri_c1_concentration;
cell_absorption = zeros(1,init_size);

nutrient_history{1} = nutrient;
absorb_history{1} = cell_absorption;

for t = 1:time
    
    % Length of colony at t
    len_fix = size(colony,2);
    
% ======================================================================
% Current nutrient allocation
% ======================================================================

    var_nutri_size = length(nutrient);

    nutrient = nutrient_distribution(var_nutri_size, ...
        nutri_movement_prob, nutrient);
    

    for i = 1:var_nutri_size
        if nutrient(i) ~= 0
            x = rand(1,nutrient(i));
            absorb_val = sum(lt(x,absorb_prob));
            cell_absorption(i) = cell_absorption(i) + absorb_val;
            nutrient(i) = nutrient(i) - absorb_val;
        end
    end
    
    x = rand(1,var_nutri_size);
    disp_temp = lt(x,replicate_prob*cell_absorption);
    displacement = cumsum(disp_temp);
    
    cell_absorption(disp_temp) = 0;
    
%    subplot(2,1,1);
%    plot(nutrient);
%    drawnow;
    
%    subplot(2,1,2);
%    plot(cell_absorption);
%    drawnow;
% Size of colony at t
    size_colony = len_fix + displacement(end);
    
% Initialise colony at t
    temp_colony = zeros(1,size_colony);
    temp_nutri = zeros(1,size_colony);
    temp_absorb = zeros(1,size_colony);
    

% ======================================================================
% Mother cell position from t-1 to t
% ======================================================================

% Old position of mother cells in colony at t-1
    position = 1:len_fix;
    
% New position of mother cells in colony at t
    position = displacement + position;
    
% Place mother cells in colony at t
    temp_colony(position) = colony;
    temp_absorb(position) = cell_absorption;
    temp_nutri(position) = nutrient;

% ======================================================================
% Daughter cell position in t
% ======================================================================
    
% Identify position of daughter cells in colony
    y = eq(temp_colony,0);
    
% Identify values of daughter cells created
    next_values = next_cell:(next_cell+displacement(end)-1);
    
% Place daughter cells in colony at t
    temp_colony(y) = next_values;
    
% Lowest cell number yet to be used
    next_cell = next_cell + displacement(end);   
    
% ======================================================================
% Length and colony history storage
% ======================================================================   
    
    colony = temp_colony;
    
    colony_history{t+1} = colony;
    len = [len, numel(colony)];
    
%    plot(len);
%    drawnow;
    
    nutrient = temp_nutri;
    nutrient(1) = nutri_c1_concentration;
    nutrient_history{t+1} = nutrient;
    
    cell_absorption = temp_absorb;
    absorb_history{t+1} = cell_absorption;
    
% ======================================================================
% Generating pathline plot
% ======================================================================     

    if pathline
        for j = 1:size(pathline_vals,2)
            
% Do not include position of values not in colony at t
            if pathline_vals(j)>size(colony,2)
                path(j,t+1) = NaN;
             
% Find position of values requested at t
            else
                path(j,t+1) = find(colony == pathline_vals(j));
            end
        end

        final_path = path;
% Pathline plot
%        plot(path,0:time,'-o')
    end
end