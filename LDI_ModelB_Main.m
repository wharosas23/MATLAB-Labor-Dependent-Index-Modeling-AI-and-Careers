%% LDI Model B Main Script
% Model B: Demand expansion + displacement
% LDI_{i,s} = (1 + beta_i * a_s)^(eta_i) * (1 - phi * AEI_i * a_s)
%
% Baseline comparison: LDI_baseline = 1 (no AI)

clear; clc; close all;

%% -------------------- CAREERS & AEI --------------------
careerNames = ["Mechanic", "Software Developer", "Graphic Designer"];

% Automation Exposure Index (AEI) — given/exogenous
AEI = [0.17, 0.36, 0.60];                  % 1 x C

%% -------------------- ADOPTION SCENARIOS --------------------
adoptionNames = ["Conservative", "Moderate", "Aggressive"];
a = [0.30, 0.50, 0.70];                    % 1 x S

%% -------------------- PARAMETERS (Model B) --------------------
% beta_i: productivity boost at full adoption (a = 1), career-specific
beta = [0.08, 0.30, 0.18];                 % 1 x C

% eta_i: demand expansion elasticity (how productivity translates to output demand)
eta  = [0.80, 1.30, 1.10];                 % 1 x C

% phi: displacement intensity scaling (shared)
phi  = 0.90;

%% -------------------- DIMENSIONS & SHAPES --------------------
C = numel(careerNames);
S = numel(a);

AEI_col  = AEI(:);                          % Cx1
beta_col = beta(:);                         % Cx1
eta_col  = eta(:);                          % Cx1

%% -------------------- SAFETY CHECK --------------------
phi_max = 1 / (max(AEI) * max(a));
if phi > phi_max
    error("phi=%.3f too large. Must be <= %.3f to keep (1-phi*AEI*a) nonnegative.", phi, phi_max);
end

%% -------------------- COMPUTE MODEL B LDI --------------------
% Productivity term
P = 1 + beta_col * a;% CxS
P = P.';

% Demand expansion term (NEW)
Q = P .^ (eta_col * ones(1,S)); % CxS
Q = Q.';

% Displacement term
D = 1 - phi * (AEI_col * a);                % CxS
D = D.';

% Labor Demand Index relative to baseline (no AI)
LDI = Q .* D; % CxS
LDI = LDI.';
assert(size(LDI,1) == numel(careerNames), ...
    'LDI rows must match number of careers');

%% -------------------- RESULTS TABLE --------------------
T = array2table(LDI, "VariableNames", matlab.lang.makeValidName(adoptionNames));
T.Career = careerNames(:);
T.AEI    = AEI_col;
T.beta   = beta_col;
T.eta    = eta_col;
T = movevars(T, ["Career","AEI","beta","eta"], "Before", 1);

disp("=== LDI Results (Model B; baseline = 1) ===");
disp(T);

writetable(T, "LDI_results_modelB.csv");

%% -------------------- PLOTS: Baseline + Scenarios --------------------
figure("Name","LDI by Scenario (Model B)","Color","w");
tiledlayout(1, C, "Padding","compact", "TileSpacing","compact");

for i = 1:C
    nexttile;

    y = [1, LDI(i,:)];                      % baseline + scenarios
    xlabels = ["Baseline", adoptionNames];

    bar(y);
    yline(1, "--", "Baseline (LDI=1)");

    title(careerNames(i));
    xlabel("Scenario");
    ylabel("Labor Demand Index (LDI)");
    xticklabels(xlabels);
    xtickangle(25);
    grid on;

    ymin = max(0, min(y) - 0.10);
    ymax = max(y) + 0.10;
    ylim([ymin, ymax]);
end

%% -------------------- QUICK DIAGNOSTICS --------------------
disp("                  === Diagnostics ===");
% --- ORIENTATION-SAFE SUMMARY TABLE ---
C = numel(careerNames);

if size(LDI,1) == C
    % LDI is C x S (careers x scenarios)
    meanLDI = mean(LDI,2);
    minLDI  = min(LDI,[],2);
    maxLDI  = max(LDI,[],2);

elseif size(LDI,2) == C
    % LDI is S x C (scenarios x careers)
    meanLDI = mean(LDI,1).';
    minLDI  = min(LDI,[],1).';
    maxLDI  = max(LDI,[],1).';

else
    error('LDI dimensions [%d %d] do not match num careers C=%d.', ...
          size(LDI,1), size(LDI,2), C);
end

T = table(careerNames(:), meanLDI, minLDI, maxLDI, ...
    'VariableNames', {'Career','Mean_LDI','Min_LDI','Max_LDI'});

disp(T)

%% -------------------- SENSITIVITY: eta sweep --------------------
% Requires: run_eta_sensitivity.m in the same folder
eta_delta = 0.20;
eta_step  = 0.02;

[Teta, ~] = run_eta_sensitivity( ...
    careerNames, adoptionNames, AEI, a, beta, eta, phi, ...
    eta_delta, eta_step);

disp("      === Minimum eta needed for LDI > 1 (within sweep grid) ===");
disp(Teta);

writetable(Teta, "eta_thresholds_LDI_crossing.csv");
