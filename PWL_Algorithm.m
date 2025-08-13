clc
clear
close all 

% MATLAB Implementation of the paper "A New Method for Improving the 
% Precision of Image Quality Assessment Metrics: Piecewise Linearization 
% of the Relationship Between the Metrics and Mean Opinion Scores.
%
% 25.06.2025
% C. Müge Bilsay

formatSpec = '%f'; 
currentDir = cd;
%% Read MOS and metric score values
mosdir = '\MetricScores_Datasets\metrics_values_TID2013\mos.txt';
msdir  = '\MetricScores_Datasets\metrics_values_TID2013\ssim.txt'; % define selected metric .txt file
mospath = append(currentDir, mosdir);
mspath  = append(currentDir, msdir);

metricName = 'SSIM'; % to use in plots

mosID = fopen(mospath,'r');
metricID = fopen(mspath,'r');

mosScore = fscanf(mosID, formatSpec);
metricScore = fscanf(metricID, formatSpec); 

N = numel(mosScore); % number of data pts. 

% For FR metrics typically [r_min, r_max] = [0,1]
% MOS range -> [0 8] 0 indicating the worst quality & 8 indicating the best
% quality
r_max = round(max(metricScore)); % (Eq.4)
r_min = round(min(metricScore)); % (Eq.4)
msRange = r_max-r_min; % metric score range

d_max = 8; % maxMOS
d_min = 0; % minMOS
MOSRange = d_max-d_min; 


Q = 10; % number of vertical divisions
Window.xrightIdx = r_max;
Window.width  = (r_max-r_min)/Q; 
Window.xleftIdx = Window.xrightIdx - Window.width; 

LF1_score = 1 - sqrt(1-metricScore); % For comparison purposes only, See [1] 
% LF2_score = 1 - sqrt(1-metricScore^2);


% Verical and horizontal division here. (Horizontal division first)

hIdx = 1; % temporarilty store horizontal vals.
vIdx = 1; % temporarilty store vertical vals.

S = zeros(N, Q); % S = [s_1, s_2, ... S_Q], (Eq.5)
M = zeros(N, Q); % M = [m_1, m_2, ... M_Q], (Eq.7)
R = zeros(Q,1); % number of horizontal windows for each vertical division (Eq. 9)
guardInterval = zeros(Q,1); % guard interval for each individual window. (Eq. 23)
Z_j    = [];
Zmos_j = [];

msArray  = [];
mosArray = [];
for i = 1:Q
    for ns = 1:N
        % Temporarilty obtain the particular MOS and metric score values if they are
        % within the selected window.(To be linearized later on)
        if (metricScore(ns) <= Window.xrightIdx && metricScore(ns)>= Window.xleftIdx)
            S(hIdx,i) = metricScore(ns); 
            M(hIdx,i) = mosScore(ns);
            hIdx = hIdx + 1; 
        end
    end

    % Update window positions
    Window.xleftIdx  = Window.xleftIdx  - Window.width;
    Window.xrightIdx = Window.xrightIdx - Window.width;

    % Re-set the horizontal index to 1.
    hIdx = 1;

    % For each individual window, obtain max-min values of both MOS and
    % metric score. (Eq.12)
    s_i_max = max(S(:,i)); % max metric score for a single vertical div. 
    s_i_min = min(S(:,i)); % min metric score for a single vertical div.
    M_i = M(:,i);
    m_i_max = max(M_i(M_i~=0)); % max MOS for a single vertical div. 
    m_i_min = min(M_i(M_i~=0)); % min MOS for a single vertical div.


    % WILL BE REMOVED!!!!!
    stdev = std(S(S(:,i)~=0,i)); % S(:,i) contains zero entries which may 
    % effect the std computation, to avoid that remove zeros and then
    % compute the std. 

    % Determine the number of horizontal divisions for each vertical window
    % i.e. R_1, R_2, ..., R_Q, See Eq.9
    alpha = 3; % scaling parameter
    if ~isempty(m_i_min)
        vSlope = (m_i_max-m_i_min)*alpha;
        R(i) = ceil(abs(vSlope)); %(Eq.9)
    else
        R(i) = 0;
    end

    beta = 0.04; % Eq.23
    guardInterval(i) = beta * R(i);

    % Perform division (Eqs. 18-20)
    Window.ytopIdx    = m_i_max; % Eq. 18
    Window.height     = (m_i_max - m_i_min)/R(i); % Eq.20 
    Window.ybottomIdx = Window.ytopIdx - Window.height; % Eq.19

    % Perform horizontal division
    nonzeroEntries = S(S(:,i)~=0,i);
    Z_j    = [Z_j, zeros(N,R(i))]; %#ok<AGROW>
    Zmos_j = [Zmos_j, zeros(N,R(i))]; %#ok<AGROW>
    for j = 1:R(i)
        for ns = 1:numel(S(:,i))
            % Set Z_idx
            if i == 1, Z_idx = j; else, Z_idx = sum(R(1:i-1))+j; end

            % Obtain metric score and MOS values within each individual
            % window
            if M_i(ns) <= Window.ytopIdx && M_i(ns) >= Window.ybottomIdx
                Z_j(vIdx, Z_idx)    = S(ns,i);
                Zmos_j(vIdx, Z_idx) = M_i(ns);
                vIdx = vIdx + 1;
            end
        end

        % Update Window Positions
        Window.ytopIdx = Window.ytopIdx - Window.height;
        Window.ybottomIdx = Window.ybottomIdx - Window.height;

        % Re-set the vertical index to 1.
        vIdx = 1;
        
        nonzero_ms = Z_j(Z_j(:,Z_idx)~=0,Z_idx);
        nonzero_mos = Zmos_j(Zmos_j(:,Z_idx)~=0,Z_idx);

        % Fit a linear function to each window using least squares.
        n = numel(nonzero_ms);
        a0 = (n*sum(nonzero_ms.*nonzero_mos) - sum(nonzero_ms)*sum(nonzero_mos)) ...
            / (n*sum(nonzero_ms.^2) - (sum(nonzero_ms))^2);
        a1 = (sum(nonzero_ms)-sum(nonzero_mos)*a0)/n;

        minVal = min(nonzero_ms); 
        maxVal = max(nonzero_ms);

        d_min_ij = (Window.ybottomIdx-guardInterval(i))*msRange/MOSRange; % Eq. 21
        d_max_ij = (Window.ytopIdx   +guardInterval(i))*msRange/MOSRange; % Eq. 22

        if (Window.ytopIdx+guardInterval) > d_max
            d_max_ij = d_max*msRange/MOSRange; % Upper limit is exceeded
        elseif (Window.ybottomIdx-guardInterval) < d_min
            d_min_ij = d_min*msRange/MOSRange; % Lower limit is exceeded
        end

        % Perform normalization
        if numel(nonzero_ms) > 1
            normalized = (nonzero_ms-minVal)/(maxVal-minVal)*(d_max_ij-d_min_ij) + d_min_ij; % Eq. 27
        else
            normalized = nonzero_ms;
        end

        if any(isnan(normalized))
            normalized = [];
            nonzero_ms = [];
            nonzero_mos = [];
        end
        msArray  = [msArray ; normalized]; %#ok<AGROW> % Store all the normalized metric score values
        mosArray = [mosArray; nonzero_mos]; %#ok<AGROW> % Store all the MOS values without changing them. 

        %% 5-parameter logistic function fit
        xgrid = linspace(0, r_max, N);
        yf = fiveParameterLogisticFunctionFit(metricScore, mosScore);

        %% Display window partitioning
        figure(15);
        plot(nonzero_ms, nonzero_mos, '+', ...
            'Color', 'b');
        hold on;
        if length(nonzero_ms)~= 0 %#ok<ISMT>
            rectangle('Position', [Window.xleftIdx+Window.width,  Window.ybottomIdx+Window.height, Window.width, Window.height], ...
                'EdgeColor', 'k', ...
                'LineWidth', 1, ...
                'LineStyle', '--');
            legend('Windows')
        end
        plot(xgrid, yf, 'r--',...
            'LineWidth', 2)
        legend('Images in TID2013', 'Fitted Curve', ...
            'FontSize', 18, ...
            'Location', 'northwest')
        ax = gca;
        ax.FontSize = 18;
        yticks([0 1 2 3 4 5 6 7 8]);

        xlabel(metricName, 'FontSize', 18);
        ylabel('MOS', 'FontSize',18);
        xlim([0 r_max])
        ylim([0 d_max])
        grid on; grid minor;

    end
end


%% Before PWL

figure(11)
scatter(metricScore, mosScore, 'Marker', '+')
hold on
xlabel(metricName, 'FontSize', 16)
ylabel('MOS', 'FontSize',16)
xlim([0 r_max]);
ylim([0 d_max]);
% xticks([0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
% yticks([0 1 2 3 4 5 6 7 8]);
legend('Images in TID2013','Fitted Curve',...
    'FontSize', 16,...
    'Location','northwest');
title('Before PWL')
%% After PWL
% % xgrid = linspace(0, r_max, N);
% % yf = fiveParameterLogisticFunctionFit(msArray, mosArray);

figure(12)
scatter(msArray, mosArray, ...
    'Marker', '+')
hold on;
xlabel(metricName,'FontSize', 16);
ylabel('MOS','FontSize', 16);
grid on; grid minor;
legend('Images in TID2013',...
    'FontSize', 16,...
    'Location','northwest');
title('After PWL')
