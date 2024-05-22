% == NMC ABSORB ==
%   Multiple simulations of NMC nutrient based cell replication.
%
% = INPUTS =
%   o init_size: initial size of the yeast colony
%   o time: number of discrete time periods for simulation
%   o rep: number of repetitions
%   o pathline: boolean indicating whether to graph pathline plots
%   o pathline_vals: values of cells to create pathline plots for
%   o c1_concentration: concentration of nutrient in first cell
%   o move_prob: chance of nutrient movement
%   o absorb_prob: chance of nutrient absorption by yeast
%   o replicate_prob: chance of proliferation per nutrient
%
% = OUTPUTS =
%   o multi_len: cell array containing the length of the colony for each 
%       time period, including when t = 0, for each BCS simulation.
%   o multi_col: cell array containing position of each yeast cell in the
%       colony for each time period, including when t = 0, for each BCS 
%       simulation.
%   o multi_nutrient: cell array containing history of nutrient count
%   o multi_absorb: cell array containing history of nutrient absorption
%   o Plot of cell count vs time for each simulation.
%   o Plot of average of cell count vs time.
%   o Pathline plots for indicated points in colony.
%   o Mean, variance, standard deviation, range of simulations.
%
% See also: basic_colony_size

clear

hold off
clf reset
hold on

% ======================================================================
% Initialising variables
% ======================================================================
% Variables used in basic_colony_size
init_size = 5;               % Integer
time = 350;                  % Integer
pathline = false;            % Boolean
pathline_vals = 1;           % Integer Array
c1_concentration = 50;
move_prob = 0.5;
replicate_prob = 0.05;
rep = 25;

% Concentration variable
absorb_prob = linspace(0,1,101);
abs_rep = length(absorb_prob);

% Storage for data from simulations
multi_col = cell(abs_rep,1);
multi_len = cell(abs_rep,1);
multi_nutrient = cell(abs_rep,1);
multi_absorb = cell(abs_rep,1);
mean_path = zeros(size(pathline_vals,2),time+1);

% ======================================================================
% Basic colony size loop
% ======================================================================

% Store length and colony data, and plot number of cells
for m = 1:abs_rep
    for k = 1:rep
        [length, colony_size, nutrient_history, absorb_history, final_path] = nutrient_model_colony(...
            init_size, time, c1_concentration, move_prob, absorb_prob(m), ...
            replicate_prob, pathline, pathline_vals);
                        
        multi_len{m}{k} = length;
        multi_col{m}{k} = colony_size;
        multi_nutrient{m}{k} = nutrient_history;
        multi_absorb{m}{k} = absorb_history;
    end
end

lobf = cell(2,abs_rep);
final_len = zeros(1,abs_rep);
max_curve = cell(2,abs_rep);

for i = 1:abs_rep
    mean_len = cell2mat(multi_len{i}');
    mean_len = sum(mean_len);
    mean_len = mean_len/rep;
    final_len(1,i) = mean_len(time);
%    plot(0:time, mean_len, '-', 'LineWidth', 2);
    
    [cubint, gof_int] = fit((0:time)',mean_len','poly3');
    cof = coeffvalues(cubint);
    syms x
    cubf = cof(1)*x^3+cof(2)*x^2+cof(3)*x+cof(4);
    cub1 = diff(cubf);
    cub2 = diff(cub1);
    lobf {1,i} = cub1;
    lobf {2,i} = cub2;
end

f = figure;
plot(absorb_prob,final_len,'-');

xlabel('Absorption Rate $\alpha$','Interpreter','latex')
ylabel('Colony Height')
title('Colony Growth With Variable Absorption Rate')

exportgraphics(f,'abslength.png','Resolution',500)

for i = 1:2
    for j = 1:abs_rep
        max_curve{i,j} = matlabFunction(lobf{i,j});
        max_curve{i,j} = max_curve{i,j}(1:time);
        max_curve{i,j} = abs(max_curve{i,j});
        max_curve{i,j} = max(max_curve{i,j});
    end
end

g = figure;
num = cell2mat(max_curve);
plot(absorb_prob,num(2,:),'-');

xlabel('Absorption Rate $\alpha$','Interpreter','latex')
ylabel('Maximum Curvature')
title('Maximum Curvature With Variable Absorption Rate')

exportgraphics(g,'abscurve.png','Resolution',500)

%{
xlabel('Time Steps')
ylabel('Yeast Cells')
title('Colony Growth With Variable Absorption')
legend('Absorb = 0.1','Absorb = 0.2','Absorb = 0.3','Absorb = 0.4','Absorb = 0.5')
yl = ylim;
hold off
%}