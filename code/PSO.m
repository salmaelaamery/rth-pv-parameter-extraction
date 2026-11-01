function [Best_fit, Best_pos, Convergence_curve] = PSO(Npop, Tmax, lb, ub, dim, fobj)
% =========================================================================
% PSO.m  —  Particle Swarm Optimization
% =========================================================================
% Kennedy & Eberhart (1995), avec facteur d'inertie decroissant.
%
% Interface standard :
%   [Best_fit, Best_pos, Convergence_curve] = PSO(Npop, Tmax, lb, ub, dim, fobj)
%
% Parametres PSO :
%   w_max = 0.9,  w_min = 0.4   (inertie decroissante lineairement)
%   c1 = 2.0  (cognitif),  c2 = 2.0  (social)
%   V_max = 0.2 * (ub - lb)
% =========================================================================

% --- Scalaire -> vecteur des bornes --------------------------------------
if isscalar(lb), lb = lb * ones(1, dim); end
if isscalar(ub), ub = ub * ones(1, dim); end

% --- Hyperparametres -----------------------------------------------------
w_max = 0.9;
w_min = 0.4;
c1    = 2.0;
c2    = 2.0;
Vmax  = 0.2 * (ub - lb);
Vmin  = -Vmax;

% --- Initialisation population -------------------------------------------
X = rand(Npop, dim) .* repmat(ub - lb, Npop, 1) + repmat(lb, Npop, 1);
V = zeros(Npop, dim);

% Evaluation initiale
fit = zeros(Npop, 1);
for i = 1:Npop
    fit(i) = fobj(X(i,:));
end

% Meilleurs personnels
pBest     = X;
pBest_fit = fit;

% Meilleur global
[Best_fit, idx] = min(pBest_fit);
Best_pos        = pBest(idx, :);

Convergence_curve = zeros(1, Tmax);

% --- Boucle principale ---------------------------------------------------
for t = 1:Tmax
    % Inertie decroissante lineairement
    w = w_max - (w_max - w_min) * t / Tmax;

    for i = 1:Npop
        r1 = rand(1, dim);
        r2 = rand(1, dim);

        % Mise a jour vitesse
        V(i,:) = w * V(i,:) ...
               + c1 * r1 .* (pBest(i,:) - X(i,:)) ...
               + c2 * r2 .* (Best_pos   - X(i,:));

        % Clampage vitesse
        V(i,:) = max(Vmin, min(Vmax, V(i,:)));

        % Mise a jour position
        X(i,:) = X(i,:) + V(i,:);

        % Clampage position
        X(i,:) = max(lb, min(ub, X(i,:)));

        % Evaluation
        f = fobj(X(i,:));

        % Mise a jour meilleur personnel
        if f < pBest_fit(i)
            pBest_fit(i) = f;
            pBest(i,:)   = X(i,:);
        end

        % Mise a jour meilleur global
        if f < Best_fit
            Best_fit = f;
            Best_pos = X(i,:);
        end
    end

    Convergence_curve(t) = Best_fit;
end

end
