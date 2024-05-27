% == BCS REPETITIONS ==
%   Multiple simulations of BCS basic cell replication.
%
% = INPUTS =
%   o init_size: initial size of the yeast colony
%   o time: number of discrete time periods for simulation
%   o prob: probability function for yeast cell proliferation
%   o rep: number of repetitions
%   o pathline: boolean indicating whether to graph pathline plots
%   o pathline_vals: values of cells to create pathline plots for
%
% = OUTPUTS =
%   o multi_len: cell array containing the length of the colony for each 
%       time period, including when t = 0, for each BCS simulation.
%   o multi_col: cell array containing position of each yeast cell in the
%       colony for each time period, including when t = 0, for each BCS 
%       simulation.
%   o Plot of cell count vs time for each simulation.
%   o Plot of average of cell count vs time.
%   o Pathline plots for indicated points in colony.
%   o Mean, variance, standard deviation, range of simulations.
%
% See also: basic_colony_size


clf reset
hold on

% ======================================================================
% Initialising variables
% ======================================================================
% Variables used in basic_colony_size
init_size = 5;              % Integer
time = 60;                  % Integer
rep = 100;                  % Integer
pathline_vals = 1:5;        % Integer Array

% Nutrient distribution function
% prob = @(x) (0.15);
 prob = @(x) (max(-0.025*x+0.4,0));
% prob = @(x) (    );

% Storage for data from simulations
multi_col = cell(rep,1);
multi_len = cell(rep,1);
mean_path = zeros(size(pathline_vals,2),time+1);

% ======================================================================
% Basic colony size loop
% ======================================================================

% Store length and colony data, and plot number of cells
for k = 1:rep
    [length, colony_size, final_path] = basic_colony_size(init_size, time, prob, ...
                            false, pathline_vals);
    multi_len{k} = length;
    multi_col{k} = colony_size;
    
% Determine the mean value of position for pathline plot
    mean_path = mean_path + final_path;
end


% ======================================================================
% Mean Colony Length Plot
% ======================================================================

mcl = subplot(1,2,1);
mean_len = cell2mat(multi_len);
mean_len = sum(mean_len);
mean_len = mean_len/rep;
plot(0:time, mean_len, '-k', 'LineWidth', 2);
hold on
for i = 1:rep
    plot(0:time,multi_len{i,1},'Color',[0 0 1 0.05])
    hold on
end

xlabel('Time Steps')
ylabel('Colony Height')
title('Average Colony Growth over Time')
legend('Mean Colony Height')

% ======================================================================
% Pathline Plot
% ======================================================================

mean_path = mean_path/rep;
plp = subplot(1,2,2);
plot(mean_path, 0:time, '-');
xlabel('Colony Position')
ylabel('Time Steps')
title('Pathline Plots')
legend(cellstr(num2str(pathline_vals', 'Cell %-d')),'Location','eastoutside')

% ======================================================================
% Metrics of final length of colony
% ======================================================================

for k = 1:rep
    final_len(k) = multi_len{k}(end);
end

len_average = mean(final_len);
len_var = var(final_len);
len_sd = std(final_len);
len_range = max(final_len) - min(final_len);

metrics = 'Mean = ' + string(len_average) + ', Variance = ' ...
    + string(len_var) + ', Standard Deviation = ' + string(len_sd) ...
    + ', Range = ' + string(len_range);
disp(metrics)