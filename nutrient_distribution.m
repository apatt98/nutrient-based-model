function nutri = nutrient_distribution(current_size, movement_prob, nutrient)

% == NUTRIENT DISTRIBUTION ==
%   Model of nutrient diffusion for nutrient model colony.

% = INPUTS =
%   o current_size: size of the yeast colony
%   o move_prob: chance of nutrient movement
%   o nutrient: current nutrient distribution

% = OUTPUTS =
%   o nutri: new nutrient distribution
    
% See also: nmc_repetitions, nutrient_model_colony


% ======================================================================
% Initial conditions of simulation
% ======================================================================

nutri = nutrient;

% ======================================================================
% Test nutrient model loop
% ======================================================================

% Create storage that accounts for overflow of nutrient
n_temp = [0 nutri 0];
% Create storage for marginal displacement of nutrient
n_amount = zeros(1,current_size+2);
    
% ==================================================================
% Nutrient displacement amount identification
% ==================================================================
    
for j = 1:current_size
        
    % Only checks for movement if position has nutrient
    if n_temp(j+1) ~= 0
        for k = 1:n_temp(j+1)
                
            % Generates two random values - one for movement, one for
            % direction
            x = rand(1,2);
            if x(1) < movement_prob
                    
                % Leftward movement
                if x(2) < 0.5
                    n_amount(j) = n_amount(j) + 1;
                    n_amount(j+1) = n_amount(j+1) - 1;
                   
                % Rightward movement
                else
                    n_amount(j+2) = n_amount(j+2) + 1;
                    n_amount(j+1) = n_amount(j+1) - 1;
                end
            end
                
        end
    end
        
end
    
        
% ==================================================================
% Nutrient displacement
% ==================================================================
n_temp = n_temp + n_amount;
    
% Moves nutrient back into simluation if it has moved beyond boundaries
if n_temp(1) > 0
    n_temp(2) = n_temp(2) + n_temp(1);
end
    
% As above. Enforces von Neumann boundary condition
if n_temp(current_size+2) > 0
    n_temp(current_size+1) = n_temp(current_size+1) + n_temp(current_size+2);
end
    
% Removes storage cells on either end of array
nutri = n_temp(2:current_size+1);