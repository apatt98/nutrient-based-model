% == NMC REPETITIONS ==
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
%clf reset

%subplot(3,3,1);
hold on

% ======================================================================
% Initialising variables
% ======================================================================
% Variables used in basic_colony_size
init_size = 50;               % Integer
time = 100;                  % Integer
rep = 100;                   % Integer
pathline = true;             % Boolean
pathline_vals = [1 25 50];         % Integer Array
c1_concentration = 50;
move_prob = 0.5;
absorb_prob = 0.1;
replicate_prob = 0.03;

% Storage for data from simulations
multi_col = cell(rep,1);
multi_len = cell(rep,1);
multi_nutrient = cell(rep,1);
multi_absorb = cell(rep,1);
mean_path = zeros(size(pathline_vals,2),time+1);

% ======================================================================
% Basic colony size loop
% ======================================================================

% Store length and colony data, and plot number of cells
for k = 1:rep
    [length, colony_size, nutrient_history, absorb_history, final_path] = nutrient_model_colony(...
        init_size, time, c1_concentration, move_prob, absorb_prob, ...
        replicate_prob, pathline, pathline_vals);
                        
    multi_len{k} = length;
    multi_col{k} = colony_size;
    multi_nutrient {k} = nutrient_history;
    multi_absorb {k} = absorb_history;
%    cgr = plot(0:time, length, '-', 'Color',[0 0 1 0.05]);
    
% Determine the mean value of position for pathline plot
    if pathline
        mean_path = mean_path + final_path;
    end
end

xlabel('Time Steps')
ylabel('Yeast Cells')
title('Colony Growth Repetitions')
yl = ylim;
%hold off

% ======================================================================
% Mean Colony Length Plot
% ======================================================================

hold on
%acg = subplot(3,3,2);
mean_len = cell2mat(multi_len);
mean_len = sum(mean_len);
mean_len = mean_len/rep;
plot(0:time, mean_len);
xlabel('Time Steps')
ylabel('Colony Height')
title('Colony Growth Repetitions')
plot(0,0,'-','Color',[1 1 1]);


% ======================================================================
% Pathline Plot
% ======================================================================

%{
f = figure;

if pathline
%    plp = subplot(3,3,[4 5]);
    mean_path = mean_path/rep;
    plot(mean_path, 0:time, '-');
    xlabel('Colony Position')
    ylabel('Time Steps')
    title('Pathline Plots')
    legend(cellstr(num2str(pathline_vals', 'Cell %-d')),'Location','eastoutside')
end
exportgraphics(f,'pathline.png','Resolution',500);

%}

% ======================================================================
% Metrics of final length of colony
% ======================================================================
final_len = zeros(1,rep);

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

% ======================================================================
% Cubic interpolation
% ======================================================================

%subplot(3,3,7);
[cubint, gof_int] = fit((0:time)',mean_len','poly3');
%cub = plot(cubint,0:time,mean_len);

%subplot(3,3,8);
cof = coeffvalues(cubint);
syms x
cubf = cof(1)*x^3+cof(2)*x^2+cof(3)*x+cof(4);
cub1 = diff(cubf);
cub2 = diff(cub1);
%fplot(cubf,[0 time]);
%fplot(cub1,[0 time]);
%hold on
%fplot(cub2,[0 time]);



% ======================================================================
% Exponantial interpolation
% ======================================================================

%{
subplot(3,3,8);
[expint, gof_exp] = fit((0:time)',mean_len','exp1')
ein = plot(expint,0:time,mean_len);
%}

% ======================================================================
% Cubic spline interpolation
% ======================================================================

%{
subplot(3,3,6);
xx = 0:0.25:time;
yy = spline(0:time, mean_len);
plot(0:time,mean_len,'x',xx,ppval(yy,xx))

subplot(3,3,9);
[breaks,coefs,l,n,d] = unmkpp(yy);
yy1 = mkpp(breaks,repmat(n-1:-1:1,d*l,1).*coefs(:,1:n-1),d);
plot(xx,ppval(yy1,xx),'-r')

hold on
[breaks,coefs,l,n,d] = unmkpp(yy1);
yy2 = mkpp(breaks,repmat(n-1:-1:1,d*l,1).*coefs(:,1:n-1),d);
plot(xx,ppval(yy2,xx),'-g')

max_curve = max(ppval(yy2,xx));

metrics = 'Mean = ' + string(len_average) + ', Variance = ' ...
    + string(len_var) + ', Standard Deviation = ' + string(len_sd) ...
    + ', Range = ' + string(len_range) + ', curve = ' + string(max_curve);
disp(metrics)
%}

% ======================================================================
% Individual nutrient heat map
% ======================================================================

%{
hm_size = 35;

nhm_data = NaN(1,hm_size,time+1);
for b = 1:time+1
    if size(nutrient_history{b},2) < hm_size
        nhm_data(1,1:size(nutrient_history{b},2),b) = nutrient_history{b};

    else
        temp_nhm = nutrient_history{b};
        nhm_data(1,1:hm_size,b) = temp_nhm(1:hm_size);
    end
end

fig = uifigure();
uip = uipanel(fig,'Position', [20 100 500 300]);
heatObj = heatmap(uip,nhm_data(:,:,1));
title(heatObj, 'Time = 0');
colormap(heatObj,parula);
xlabel("Cell");
n = size(nhm_data,3);
uis = uislider(fig,'Position',[50 50 450 3],...
    'Value',1,...
    'Limits',[1,n],...
    'MajorTicks',[1, 25:25:n],...
    'MinorTicks',[]);
caxis manual
caxis([0 c1_concentration])
uis.ValueChangingFcn = {@sliderChangingFcn, nhm_data, heatObj};


function sliderChangingFcn(~,event,nhm_data,heatObj)
value = round(event.Value);
heatObj.ColorData = nhm_data(:,:,value);
heatObj.Title = sprintf('Time = %d', value-1);
end
%}

% ======================================================================
% Mean nutrient heat map
% ======================================================================

%{
hm_size = 120;
mean_nhm = cell(rep,1);

for a = 1:rep
    nhm_data = NaN(1,hm_size,time+1);

    for b = 1:time+1

        if size(multi_nutrient{a}{b},2) < hm_size
            nhm_data(1,1:size(multi_nutrient{a}{b},2),b) = multi_nutrient{a}{b};
        else
            temp_nhm = multi_nutrient{a}{b};
            nhm_data(1,1:hm_size,b) = temp_nhm(1:hm_size);
        end
    end
    mean_nhm{a} = nhm_data;
end

mean_nhm = cell2mat(mean_nhm);
mean_nhm_final = nansum(mean_nhm,1)/rep;
for a = 1:hm_size
    for b = 1:time+1
        if isnan(mean_nhm(:,a,b))
            mean_nhm_final(1,a,b) = NaN;
        end
    end
end

hold off
heatObj = heatmap(mean_nhm_final(:,:,1),'xlabel','Yeast Cells',...
    'ylabel','Nutrient','title','Mean Diffusion of Nutrient over Time (120 Cells)');
colormap(heatObj,parula);
ax = gca;
ax.XDisplayLabels = nan(size(ax.XDisplayData));
ax.YDisplayLabels = nan(size(ax.YDisplayData));
n = size(mean_nhm_final,3);
caxis([0 c1_concentration])

gif('nutridiff.gif','DelayTime',0.15);
for k = 2:time+1
    set(heatObj,'ColorData',mean_nhm_final(:,:,k));
    gif
end
%}