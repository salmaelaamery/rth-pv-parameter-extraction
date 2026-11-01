function [Best_fit, Best_pos, Convergence_curve] = GWO(Npop, Tmax, lb, ub, dim, fobj)
% =========================================================================
% GWO.m  —  Grey Wolf Optimizer
% =========================================================================
% Mirjalili et al. (2014), Advances in Engineering Software.
%
% Interface standard :
%   [Best_fit, Best_pos, Convergence_curve] = GWO(Npop, Tmax, lb, ub, dim, fobj)
%
% Hierarchie : Alpha (meilleur), Beta (2e), Delta (3e), Omega (reste)
% Parametre a : decroit lineairement de 2 a 0
% =========================================================================

% --- Scalaire -> vecteur des bornes --------------------------------------
if isscalar(lb), lb = lb * ones(1, dim); end
if isscalar(ub), ub = ub * ones(1, dim); end

% --- Initialisation population -------------------------------------------
X = rand(Npop, dim) .* repmat(ub - lb, Npop, 1) + repmat(lb, Npop, 1);

% Initialisation des trois meilleurs loups
Alpha_pos = zeros(1, dim);  Alpha_score = inf;
Beta_pos  = zeros(1, dim);  Beta_score  = inf;
Delta_pos = zeros(1, dim);  Delta_score = inf;

% Evaluation initiale
for i = 1:Npop
    f = fobj(X(i,:));
    if f < Alpha_score
        Delta_score = Beta_score;  Delta_pos = Beta_pos;
        Beta_score  = Alpha_score; Beta_pos  = Alpha_pos;
        Alpha_score = f;           Alpha_pos = X(i,:);
    elseif f < Beta_score
        Delta_score = Beta_score;  Delta_pos = Beta_pos;
        Beta_score  = f;           Beta_pos  = X(i,:);
    elseif f < Delta_score
        Delta_score = f;           Delta_pos = X(i,:);
    end
end

Convergence_curve = zeros(1, Tmax);

% --- Boucle principale ---------------------------------------------------
for t = 1:Tmax
    % Facteur a : decroit lineairement de 2 a 0
    a = 2 - 2 * t / Tmax;

    for i = 1:Npop
        % --- Mise a jour position par rapport aux 3 leaders ---------------
        r1 = rand(1, dim); r2 = rand(1, dim);
        A1 = 2 * a * r1 - a;
        C1 = 2 * r2;
        D_alpha = abs(C1 .* Alpha_pos - X(i,:));
        X1 = Alpha_pos - A1 .* D_alpha;

        r1 = rand(1, dim); r2 = rand(1, dim);
        A2 = 2 * a * r1 - a;
        C2 = 2 * r2;
        D_beta = abs(C2 .* Beta_pos - X(i,:));
        X2 = Beta_pos - A2 .* D_beta;

        r1 = rand(1, dim); r2 = rand(1, dim);
        A3 = 2 * a * r1 - a;
        C3 = 2 * r2;
        D_delta = abs(C3 .* Delta_pos - X(i,:));
        X3 = Delta_pos - A3 .* D_delta;

        X(i,:) = (X1 + X2 + X3) / 3;

        % Clampage
        X(i,:) = max(lb, min(ub, X(i,:)));

        % Evaluation
        f = fobj(X(i,:));

        % Mise a jour hierarchie
        if f < Alpha_score
            Delta_score = Beta_score;  Delta_pos = Beta_pos;
            Beta_score  = Alpha_score; Beta_pos  = Alpha_pos;
            Alpha_score = f;           Alpha_pos = X(i,:);
        elseif f < Beta_score
            Delta_score = Beta_score;  Delta_pos = Beta_pos;
            Beta_score  = f;           Beta_pos  = X(i,:);
        elseif f < Delta_score
            Delta_score = f;           Delta_pos = X(i,:);
        end
    end

    Convergence_curve(t) = Alpha_score;
end

Best_fit = Alpha_score;
Best_pos = Alpha_pos;

end
