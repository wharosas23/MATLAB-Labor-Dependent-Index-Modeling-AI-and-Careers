function [Teta, LDI_sweep] = run_eta_sensitivity( ...
    careerNames, adoptionNames, AEI, a, beta, eta, phi, ...
    eta_delta, eta_step)

% RUN_ETA_SENSITIVITY
% Sensitivity analysis for Model B:
% LDI_{i,s} = (1 + beta_i a_s)^(eta_i) * (1 - phi AEI_i a_s)
%
% Outputs:
%   Teta       : table of minimum eta needed for LDI>1
%   LDI_sweep  : C x S x K array of LDI values

%% --- Dimensions & reshaping ---
C = numel(careerNames);
S = numel(a);

AEI_col  = AEI(:);
beta_col = beta(:);
eta_base = eta(:);

%% --- Eta grid ---
eta_grid = -eta_delta:eta_step:eta_delta;
K = numel(eta_grid);

LDI_sweep = zeros(C, S, K);

%% --- Precompute terms ---
P = 1 + beta_col * a;        % productivity term (CxS)
D = 1 - phi * (AEI_col * a); % displacement term (CxS)

%% --- Sweep ---
for k = 1:K
    eta_k = eta_base + eta_grid(k);        % shift all etas together
    Qk = P .^ (eta_k * ones(1,S));         % demand expansion
    LDI_sweep(:,:,k) = Qk .* D;
end

%% --- Minimum eta needed for LDI > 1 ---
eta_min_needed = nan(C,S);

for i = 1:C
    for s = 1:S
        idx = find(squeeze(LDI_sweep(i,s,:)) > 1, 1, "first");
        if ~isempty(idx)
            eta_min_needed(i,s) = eta_base(i) + eta_grid(idx);
        end
    end
end

Teta = array2table(eta_min_needed, ...
    "VariableNames", matlab.lang.makeValidName(adoptionNames));
Teta.Career = careerNames(:);
Teta = movevars(Teta, "Career", "Before", 1);

end
