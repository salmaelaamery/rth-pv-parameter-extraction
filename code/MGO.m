function [Best_fit, Best_pos, Convergence_curve] = MGO(Npop, Tmax, lb, ub, dim, fobj)
% =========================================================================
% MGO.m  —  Mountain Gazelle Optimizer
% =========================================================================
% Abdollahzadeh et al. (2022), Advances in Engineering Software.
%
% Interface standard :
%   [Best_fit, Best_pos, Convergence_curve] = MGO(Npop, Tmax, lb, ub, dim, fobj)
%
% Inspire du comportement social des gazelles de montagne :
%   - Harems independants (exploration locale autour du mâle dominant)
%   - Gazelles solitaires (exploration globale)
%   - Migration (exploitation des meilleures zones)
%   - Recherche de nourriture (exploitation fine)
% =========================================================================

% --- Scalaire -> vecteur des bornes --------------------------------------
if isscalar(lb), lb = lb * ones(1, dim); end
if isscalar(ub), ub = ub * ones(1, dim); end

% --- Initialisation population -------------------------------------------
X = rand(Npop, dim) .* repmat(ub - lb, Npop, 1) + repmat(lb, Npop, 1);

fit = zeros(Npop, 1);
for i = 1:Npop
    fit(i) = fobj(X(i,:));
end

[Best_fit, idx] = min(fit);
Best_pos = X(idx, :);

Convergence_curve = zeros(1, Tmax);

% --- Fonction Levy -------------------------------------------------------
    function L = levy_mgo(n, m, beta)
        num   = gamma(1 + beta) * sin(pi * beta / 2);
        den   = gamma((1 + beta)/2) * beta * 2^((beta-1)/2);
        sigma = (num / den)^(1 / beta);
        u = randn(n, m) * sigma;
        v = randn(n, m);
        L = u ./ (abs(v).^(1 / beta));
    end

% --- Boucle principale ---------------------------------------------------
for t = 1:Tmax
    % Parametre de controle (diminue avec t)
    M = exp(-t / Tmax);

    X_mean = mean(X, 1);

    for i = 1:Npop
        X_new = zeros(1, dim);
        r1 = rand; r2 = rand; r3 = rand; r4 = rand;

        % Selection stochastique du comportement
        choice = ceil(4 * rand);

        switch choice
            case 1
                % Harem — exploration locale autour de Best_pos
                cof = 1 / (dim + 1);
                Xr  = lb + rand(1,dim) .* (ub - lb);  % position aleatoire
                X_new = Best_pos + ...
                        randn(1,dim) .* abs(cof * (Best_pos - Xr)) + ...
                        r1 * (Best_pos - X_mean);

            case 2
                % Gazelle solitaire — exploration globale avec Levy
                LF = levy_mgo(1, dim, 1.5);
                X_new = Best_pos .* LF + ...
                        lb + rand(1,dim) .* (ub - lb);

            case 3
                % Migration — deplacement vers Best_pos avec perturbation
                A = 2 * r2 - 1;   % dans [-1, 1]
                X_new = Best_pos - ...
                        abs(rand(1,dim) .* Best_pos - X(i,:)) .* ...
                        (A + M * randn(1,dim));

            case 4
                % Recherche de nourriture — exploitation fine
                Xr1 = X(randi(Npop),:);
                Xr2 = X(randi(Npop),:);
                X_new = X(i,:) + ...
                        r3 .* (Xr1 - X(i,:)) + ...
                        r4 .* (Xr2 - Best_pos);
        end

        % Clampage
        X_new = max(lb, min(ub, X_new));

        % Evaluation & mise a jour greedy
        f_new = fobj(X_new);
        if f_new < fit(i)
            X(i,:) = X_new;
            fit(i) = f_new;
        end

        if fit(i) < Best_fit
            Best_fit = fit(i);
            Best_pos = X(i,:);
        end
    end

    Convergence_curve(t) = Best_fit;
end

end
