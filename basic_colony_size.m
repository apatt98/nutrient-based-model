function [len, colony_history, final_path] = basic_colony_size(init_size, time, ...
                                 prob, pathline, pathline_vals)
% == BASIC COLONY SIZE ==
%   One simulation of basic cell replication. Inserts yeast cells into
%   a colony based on a probability function and discrete time steps.
%
% = INPUTS =
%   o init_size: initial size of the yeast colony
%   o time: number of discrete time periods for simulation
%   o prob: probability function for yeast cell proliferation
%   o pathline: boolean to create pathline graph for given values; defaults
%       to false.
%   o pathline_vals: array of yeast cells for pathline graph; defaults to
%       yeast cell 1.
%
% = OUTPUTS =
%   o length: an array containing the size of the colony for each time
%       period, including when t = 0.
%   o colony_history: cell array containing the position of each yeast
%       cell in the colony for each time period, including when t = 0. 
%   o Plot of position of specific yeast cells over time (optional).
%       
% See also: bcs_repetitions



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

for t = 1:time
    
% Length of colony at t
    len_fix = size(colony,2);
    
% ======================================================================
% Yeast cell proliferation chance
% ======================================================================

% Probability of proliferation for each cell in colony from t-1
    probability = prob(1:len_fix);
    
% Random value for each cell in colony to determine proliferation
    proliferation = rand(1,len_fix);
    
% Displacement of cells in colony from t-1 to t
    displacement = cumsum(lt(proliferation,probability));
    
% Size of colony at t
    size_colony = len_fix + displacement(end);
    
% Initialise colony at t
    temp_colony = zeros(1,size_colony);
    
    
% ======================================================================
% Mother cell position from t-1 to t
% ======================================================================

% Old position of mother cells in colony at t-1
    position = 1:len_fix;
    
% New position of mother cells in colony at t
    position = displacement + position;
    
% Place mother cells in colony at t
    temp_colony(position) = colony;

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