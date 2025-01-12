
%% Introduction
% This script implements a proposed statistical method for identifying drought and heatwave events from climate data. 
% It includes data pre-processing, calculation of SPI and SHI, events identification, and visualization.

% Please cite this paper if you use this method:  
% Shan, B., Verhoest, N. E. C., and De Baets, B.: Identification of compound drought and heatwave events on a daily scale and across four seasons, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2023-147, 2023.

% Auther: Baoying Shan, Gent University, (baoying.shan@ugent.be)
% Date: 2024.03.24


% This script includes
% part 1: Load functions and data
% part 2: Calculate SPI and SHI
% Part 3: Removal and merging processes
% Part 4: Identify drought and heatwave events
% Part 5: Identify compound drought and heatwave events
% Part 6: Plotting

%% Part 1: Load functions and data

addpath(genpath('./02_src'));
data = readtable('01_data/CGI_Weather_MATLAB_Formatted_2.csv')

% which saved Pre(the daily precipitation), DMT(daily mean temperature), and Date (corresponding year, month, day) 

% Data pre-processing:
% Convert 366 days to 365 days
mm = find(Date(:,2) == 2 & Date(:,3) == 29);
Pre = pre; Pre(mm - 1, :) = pre(mm - 1, :) + pre(mm, :);
DMT = dmt; DMT(mm - 1, :) = (dmt(mm - 1, :) + dmt(mm, :)) / 2;
Date(mm, :) = []; Pre(mm, :) = []; DMT(mm, :) = [];
t = datetime(Date);
years = Date(end, 1) - Date(1, 1) + 1;

%% Part 2: Calculate SPI and SHI

NSP = 30; % parameter for non-stationarity
scale_p = 30; % accumulation period for droughts
scale_T = 3; % accumulation period for heatwaves

% Parametric method to calculate SPI and SHI
SPI = SPI_best_dist(Date, Pre, scale_p, NSP);
[SHI] = SHI_best_dist(Date, DMT, scale_T, NSP);

% nonparametric method to calculate SPI and SHI to improve the efficiency
%SPI = SI_nonparametric(Date, pre_one_grid, scale_p, NSP);
%SHI = SI_nonparametric(Date, dmat_one_grid, scale_T, NSP);

%% Part 3: Removal and merging processes
% Parameter setting
start_th_d=-1; end_th_d=-1; % pre-identification thresholds for droughts
start_th_h=1;  end_th_h=1; %  pre-identification thresholds for heatwaves
start_th_p=1;  end_th_p=1; %  pre-identification thresholds for wet events
start_th_c=-1;  end_th_c=-1; %  pre-identification thresholds for heatwaves
p=0.05;
events_number_control=0.1; events_days_control=120;

% Get removal and merging thresholds for droughts
cc =  remo_merg("d",Date, SPI, ...
    scale_p, p, start_th_d, end_th_d, events_number_control, events_days_control); % "d" represents droughts
drought_removal_threshold=cc(1)  
drought_merging_threshold=cc(2)

% Get removal and merging thresholds for heatwaves
cc = remo_merg("h",Date, SHI, ...
    scale_T, p, start_th_h, end_th_h, events_number_control, events_days_control); % "h" represnets heatwaves
heatwave_removal_threshold=cc(1) 
heatwave_merging_threshold=cc(2)

% Get removal and merging thresholds for wet events
cc = remo_merg("p", Date, SPI, ...
    scale_p, p, start_th_p, end_th_p, events_number_control, events_days_control); % "p" represents wet events
wet_removal_threshold = cc(1);  
wet_merging_threshold = cc(2);

% Get removal and merging thresholds for cold waves
cc = remo_merg("c",Date, SHI, ...
    scale_T, p, start_th_c, end_th_c, events_number_control, events_days_control); % "c" represnets colwaves
coldwave_removal_threshold=cc(1) 
coldwave_merging_threshold=cc(2)


%% Part 4: Identify drought and heatwave events
drought_daily = PRM_extreme_identification("d", Date, SPI,start_th_d,end_th_d,...
    drought_removal_threshold, drought_merging_threshold ); % the ninth column records if this day is in a drought (1 for yes, 0 for no)
heatwave_daily = PRM_extreme_identification("h", Date, SHI, start_th_h, end_th_h,...
    heatwave_removal_threshold, heatwave_merging_threshold); % the ninth column records if this day is in a heatwave (1 for yes, 0 for no)
wet_daily = PRM_extreme_identification("p", Date, SPI,start_th_p,end_th_p,...
    wet_removal_threshold, wet_merging_threshold ); % the ninth column records if this day is in a wet event (1 for yes, 0 for no)
coldwave_daily = PRM_extreme_identification("c", Date, SHI, start_th_c, end_th_c,...
    coldwave_removal_threshold, coldwave_merging_threshold); % the ninth column records if this day is in a coldwave (1 for yes, 0 for no)

drought_events = daily_2_events(drought_daily, "d"); % drought events
heatwave_events = daily_2_events(heatwave_daily, "h"); % heatwave events
wet_events = daily_2_events(wet_daily, "p"); % wet events
coldwave_events = daily_2_events(coldwave_daily, "c"); % coldwave events

%% Part 5: Identify compound events
types=["d-and-h", "d-or-h", "d-cond-h", "h-cond-d", "d-and-c", "d-or-c", "d-cond-c", "c-cond-d", "p-and-c", "p-or-c", "p-cond-c", "c-cond-p", "d-and-c", "d-or-c", "d-cond-c", "c-cond-d", "p-and-h", "p-or-h", "p-cond-h", "h-cond-p"];
% Initialize storage for compound event results
compound_events = cell(length(types), 1); % Cell array to hold results for each type

% Loop through each type and process compound events
for i = 1:length(types)
    type = types(i);
    
    % Determine which datasets and mode to use based on the type
    switch type
        case "d-and-h" % Drought and heatwaves
            [~, compound_daily] = identify_compound(drought_daily, heatwave_daily, 1);
        case "d-or-h" % Drought or heatwaves
            [~, compound_daily] = identify_compound(drought_daily, heatwave_daily, 2);
        case "d-cond-h" % Drought conditioned on heatwaves
            [~, compound_daily] = identify_compound(drought_daily, heatwave_daily, 3);
        case "h-cond-d" % Heatwaves conditioned on drought
            [~, compound_daily] = identify_compound(drought_daily, heatwave_daily, 4);
        case "d-and-c" % Drought and coldwaves
            [~, compound_daily] = identify_compound(drought_daily, coldwave_daily, 1);
        case "d-or-c" % Drought or coldwaves
            [~, compound_daily] = identify_compound(drought_daily, coldwave_daily, 2);
        case "d-cond-c" % Drought conditioned on coldwaves
            [~, compound_daily] = identify_compound(drought_daily, coldwave_daily, 3);
        case "c-cond-d" % Coldwaves conditioned on drought
            [~, compound_daily] = identify_compound(drought_daily, coldwave_daily, 4);
        case "p-and-c" % Wet events and coldwaves
            [~, compound_daily] = identify_compound(wet_daily, coldwave_daily, 1);
        case "p-or-c" % Wet events or coldwaves
            [~, compound_daily] = identify_compound(wet_daily, coldwave_daily, 2);
        case "p-cond-c" % Wet events conditioned on coldwaves
            [~, compound_daily] = identify_compound(wet_daily, coldwave_daily, 3);
        case "c-cond-p" % Coldwaves conditioned on wet events
            [~, compound_daily] = identify_compound(wet_daily, coldwave_daily, 4);
        case "p-and-h" % Wet events and heatwaves
            [~, compound_daily] = identify_compound(wet_daily, heatwave_daily, 1);
        case "p-or-h" % Wet events or heatwaves
            [~, compound_daily] = identify_compound(wet_daily, heatwave_daily, 2);
        case "p-cond-h" % Wet events conditioned on heatwaves
            [~, compound_daily] = identify_compound(wet_daily, heatwave_daily, 3);
        case "h-cond-p" % Heatwaves conditioned on wet events
            [~, compound_daily] = identify_compound(wet_daily, heatwave_daily, 4);
        otherwise
            error('Unknown compound event type: %s', type);
    end
    
    % Store the result in the cell array
    compound_events{i} = compound_daily;
end

%% Part 6: Plotting

% Define the year range for analysis
year_range = [2021, 2021];
time_lim = [t(365 * (year_range(1) - Date(1, 1)) + 1), t(365 * (year_range(2) - Date(1, 1)) + 365)];

figure(1);

% Plot SPI and drought identification results
ax(1) = subplot(10, 1, [1:2]);
hold on;
cc = drought_daily(:, 4);
cc(drought_daily(:, end - 1) == 0) = nan;
ba = bar(t, cc, 'FaceColor', [1, 0.8, 0]);
ba(1).BarWidth = 1;
plot(t, SPI, '.-', 'color', [0.85, 0.33, 0.10], 'LineWidth', 0.8);
plot(t, start_th_d * ones(length(t), 1), 'r--', 'LineWidth', 0.5);
xlim(time_lim);
hold off;
ylabel('SPI');
grid on;
title('Drought');
datetick('x', 'yyyy-mm');
set(gca, 'xticklabel', []);
box on;

% Plot SHI and heatwave identification results
ax(2) = subplot(10, 1, [3:4]);
hold on;
cc = heatwave_daily(:, 4);
cc(heatwave_daily(:, end - 1) == 0) = nan;
ba = bar(t, cc, 'FaceColor', [0.6, 0.8, 0.5]);
ba(1).BarWidth = 1;
plot(t, SHI, '.-', 'color', [0.47, 0.67, 0.19], 'LineWidth', 0.8);
plot(t, start_th_h * ones(length(t), 1), 'r--', 'LineWidth', 0.5);
hold off;
xlim(time_lim);
ylim([-3, 3]);
title('Heatwave');
ylabel('SHI');
grid on;
datetick('x', 'yyyy-mm');
set(gca, 'xticklabel', []);
box on;

% Plot coldwave identification results
ax(3) = subplot(10, 1, [5:6]);
hold on;
cc = coldwave_daily(:, 4);
cc(coldwave_daily(:, end - 1) == 0) = nan;
ba = bar(t, cc, 'FaceColor', [0.3, 0.5, 0.9]);
ba(1).BarWidth = 1;
plot(t, SHI, '.-', 'color', [0.2, 0.4, 0.8], 'LineWidth', 0.8);
plot(t, start_th_c * ones(length(t), 1), 'r--', 'LineWidth', 0.5);
hold off;
xlim(time_lim);
ylim([-3, 3]);
title('Coldwave');
ylabel('CWI');
grid on;
datetick('x', 'yyyy-mm');
set(gca, 'xticklabel', []);
box on;

% Plot wet events identification results
ax(4) = subplot(10, 1, [7:8]);
hold on;
cc = wet_daily(:, 4);
cc(wet_daily(:, end - 1) == 0) = nan;
ba = bar(t, cc, 'FaceColor', [0.2, 0.8, 0.8]);
ba(1).BarWidth = 1;
plot(t, SPI, '.-', 'color', [0.1, 0.7, 0.7], 'LineWidth', 0.8);
plot(t, start_th_p * ones(length(t), 1), 'r--', 'LineWidth', 0.5);
hold off;
xlim(time_lim);
title('Wet Events');
ylabel('SPI');
grid on;
datetick('x', 'yyyy-mm');
set(gca, 'xticklabel', []);
box on;

% Plot compound events for all types
for i = 1:length(types)
    type = i;
    compound = identify_compound(drought_daily, heatwave_daily, type); % Update for the specific compound combination

    ax(4 + i) = subplot(10, 1, 8 + i); % For compound period
    ba = bar(t, compound(:, 1));
    ba(1).BarWidth = 2;
    ba(1).FaceColor = [0.5, 0.5, 0.5];
    xlim(time_lim);
    yticks([0, 1]);
    ylabel('Compound');
    set(gca, 'yticklabel', []);
    title(types(type));
    grid on;

    if type < length(types)
        datetick('x', 'yyyy-mm');
        set(gca, 'xticklabel', []);
    else
        datetick('x', 'yyyy-mm');
        xlabel('Date');
    end
end

linkaxes(ax, 'x');
set(gcf, 'Position', [150, 120, 1000, 800]);

% Save the plot
% print(gcf, 'identification_results.png', '-dpng', '-r600');

%% Part 6: Plotting for 2021 (April 26th to October 29th)

% Define the time range
start_date = datetime(2021, 4, 26);
end_date = datetime(2021, 10, 29);
time_lim = [start_date, end_date];

figure(1);

% Plot SPI and drought identification results
ax(1) = subplot(5, 1, 1); % 5 subplots, this is the 1st
hold on;
cc = drought_daily(:, 4);
cc(drought_daily(:, end - 1) == 0) = nan;
ba = bar(t, cc, 'FaceColor', [1, 0.8, 0]);
ba(1).BarWidth = 1;
plot(t, SPI, '.-', 'color', [0.85, 0.33, 0.10], 'LineWidth', 0.8);
plot(t, start_th_d * ones(length(t), 1), 'r--', 'LineWidth', 0.5);
xlim(time_lim);
hold off;
ylabel('SPI');
grid on;
title('Drought');
datetick('x', 'yyyy-mm');
set(gca, 'xticklabel', []);
box on;

% Plot SHI and heatwave identification results
ax(2) = subplot(5, 1, 2); % 5 subplots, this is the 2nd
hold on;
cc = heatwave_daily(:, 4);
cc(heatwave_daily(:, end - 1) == 0) = nan;
ba = bar(t, cc, 'FaceColor', [0.6, 0.8, 0.5]);
ba(1).BarWidth = 1;
plot(t, SHI, '.-', 'color', [0.47, 0.67, 0.19], 'LineWidth', 0.8);
plot(t, start_th_h * ones(length(t), 1), 'r--', 'LineWidth', 0.5);
hold off;
xlim(time_lim);
ylim([-3, 3]);
title('Heatwave');
ylabel('SHI');
grid on;
datetick('x', 'yyyy-mm');
set(gca, 'xticklabel', []);
box on;

% Plot SHI and coldwave identification results
ax(3) = subplot(5, 1, 3); % 5 subplots, this is the 3rd
hold on;
cc = coldwave_daily(:, 4);
cc(coldwave_daily(:, end - 1) == 0) = nan;
ba = bar(t, cc, 'FaceColor', [0.3, 0.5, 0.9]);
ba(1).BarWidth = 1;
plot(t, SHI, '.-', 'color', [0.2, 0.4, 0.8], 'LineWidth', 0.8);
plot(t, start_th_c * ones(length(t), 1), 'r--', 'LineWidth', 0.5);
hold off;
xlim(time_lim);
ylim([-3, 3]);
title('Coldwave');
ylabel('SHI');
grid on;
datetick('x', 'yyyy-mm');
set(gca, 'xticklabel', []);
box on;

% Plot SPI and wet events identification results
ax(4) = subplot(5, 1, 4); % 5 subplots, this is the 4th
hold on;
cc = wet_daily(:, 4);
cc(wet_daily(:, end - 1) == 0) = nan;
ba = bar(t, cc, 'FaceColor', [0.2, 0.8, 0.8]);
ba(1).BarWidth = 1;
plot(t, SPI, '.-', 'color', [0.1, 0.7, 0.7], 'LineWidth', 0.8);
plot(t, start_th_p * ones(length(t), 1), 'r--', 'LineWidth', 0.5);
hold off;
xlim(time_lim);
title('Wet Events');
ylabel('SPI');
grid on;
datetick('x', 'yyyy-mm');
set(gca, 'xticklabel', []);
box on;

% Plot p-and-c compound events
[~, compound_daily] = identify_compound(wet_daily, coldwave_daily, 1);
compound_filtered = compound_daily(t >= start_date & t <= end_date, :); % Filter for date range
dates_filtered = t(t >= start_date & t <= end_date);

ax(5) = subplot(5, 1, 5); % 5 subplots, this is the 5th
hold on;
bar(dates_filtered, compound_filtered(:, end), 'FaceColor', [0.5, 0.5, 0.5], 'BarWidth', 1);
xlim(time_lim);
ylim([-0.5, 1.5]); % Binary values (0 or 1)
yticks([0, 1]);
ylabel('P-and-C');
xlabel('Date');
title('P-and-C Compound Events');
grid on;
datetick('x', 'yyyy-mm-dd', 'keeplimits');
hold off;

% Link axes for synchronization
linkaxes(ax, 'x');

% Adjust figure size
set(gcf, 'Position', [150, 120, 1000, 800]);

% Save the plot
% print(gcf, 'compound_events_2021.png', '-dpng', '-r600');
