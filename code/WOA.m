function [Best_fit, Best_pos, Convergence_curve] = WOA(Npop, Tmax, lb, ub, dim, fobj)
% =========================================================================
% WOA.m  —  Whale Optimization Algorithm
% =========================================================================
% Mirjalili & Lewis (2016), Advances in Engineering Software.
%
% Interface standard :
%   [Best_fit, Best_pos, Convergence_curve] = WOA(Npop, Tmax, lb, ub, dim, fobj)
%
% Mecanismes :
%   - Encerclement en reduction (shrinking encircling)
%   - Attaque en spirale (bubble-net feeding)
%   - Recherche de proie aleatoire (exploration)
%   b = 1 (constante spirale), l dans [-1, 1]
% =========================================================================

% --- Scalaire -> vecteur des bornes --------------------------------------
if isscalar(lb), lb = lb * ones(1, dim); end
if isscalar(ub), ub = ub * ones(1, dim); end

% --- Initialisation population -------------------------------------------
X = rand(Npop, dim) .* repmat(ub - lb, Npop, 1) + repmat(lb, Npop, 1);

% Evaluation initiale
fit = zeros(Npop, 1);
for i = 1:Npop
    fit(i) = fobj(X(i,:));
end

[Best_fit, idx] = min(fit);
Best_pos = X(idx, :);

b = 1;   % constante spirale

Convergence_curve = zeros(1, Tmax);

% --- Boucle principale ---------------------------------------------------
for t = 1:Tmax
    a  = 2 - 2 * t / Tmax;   % decroit de 2 a 0
    a2 = -1 + (-1) * t / Tmax; % decroit de -1 a -2 (pour l)

    for i = 1:Npop
        r1 = rand;
        r2 = rand;
        A  = 2 * a * r1 - a;
        C  = 2 * r2;
        l  = (a2 - 1) * rand + 1;  % l dans [a2, 1]
        p  = rand;

        if p < 0.5
            if abs(A) < 1
                % Encerclement en reduction
                D = abs(C .* Best_pos - X(i,:));
                X(i,:) = Best_pos - A .* D;
            else
                % Recherche aleatoire
                rand_idx = randi(Npop);
                X_rand = X(rand_idx, :);
                D = abs(C .* X_rand - X(i,:));
                X(i,:) = X_rand - A .* D;
            end
        else
            % Attaque en spirale
            D_leader = abs(Best_pos - X(i,:));
            X(i,:) = D_leader .* exp(b * l) .* cos(2 * pi * l) + Best_pos;
        end

        % Clampage
        X(i,:) = max(lb, min(ub, X(i,:)));

        % Evaluation
        f = fobj(X(i,:));

        if f < Best_fit
            Best_fit = f;
            Best_pos = X(i,:);
        end
        fit(i) = f;
    end

    Convergence_curve(t) = Best_fit;
end

end
