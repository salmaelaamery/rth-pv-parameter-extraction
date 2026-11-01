function [Best_fit, Best_pos, Convergence_curve] = AVOA(Npop, Tmax, lb, ub, dim, fobj)
% =========================================================================
% AVOA.m  —  African Vultures Optimization Algorithm
% =========================================================================
% Abdollahzadeh et al. (2021), Computers & Industrial Engineering.
% doi: 10.1016/j.cie.2021.107408
%
% Interface standard :
%   [Best_fit, Best_pos, Convergence_curve] = AVOA(Npop, Tmax, lb, ub, dim, fobj)
%
% Parametres :
%   p1=0.6, p2=0.4, p3=0.6  (probabilites de selection de phase)
%   alpha=0.8, Betha=0.2     (ponderation des deux meilleurs vautours)
%
% CORRECTION v2 :
%   Formule satiation : F = (2*rand+1)*(1 - t/Tmax)
%   => |F| decroit de [1,3] vers 0 => transition exploration->exploitation
%   La formule precedente avec sin() gardait |F|>=1 (exploration permanente)
% =========================================================================

% --- Scalaire -> vecteur des bornes --------------------------------------
if isscalar(lb), lb = lb * ones(1, dim); end
if isscalar(ub), ub = ub * ones(1, dim); end

% --- Hyperparametres (Abdollahzadeh 2021) --------------------------------
p1    = 0.6;
p2    = 0.4;
p3    = 0.6;
alpha = 0.8;
Betha = 0.2;

% --- Initialisation population -------------------------------------------
X = rand(Npop, dim) .* repmat(ub - lb, Npop, 1) + repmat(lb, Npop, 1);

fit = zeros(Npop, 1);
for i = 1:Npop
    fit(i) = fobj(X(i,:));
end

% Deux meilleurs vautours (classement initial)
[sorted_fit, sorted_idx] = sort(fit);
Best1_pos = X(sorted_idx(1), :);  Best1_fit = sorted_fit(1);
Best2_pos = X(sorted_idx(2), :);  Best2_fit = sorted_fit(2);

Convergence_curve = zeros(1, Tmax);

% --- Fonction Levy -------------------------------------------------------
    function L = levy_flight(n, m, beta)
        num   = gamma(1 + beta) * sin(pi * beta / 2);
        den   = gamma((1 + beta)/2) * beta * 2^((beta-1)/2);
        sigma = (num / den)^(1 / beta);
        u = randn(n, m) * sigma;
        v = randn(n, m);
        L = u ./ (abs(v).^(1 / beta));
    end

% --- Boucle principale ---------------------------------------------------
for t = 1:Tmax

    % Energie de satiation : decroit de [1,3] vers 0
    % CORRECTION : formule lineaire sans terme sin (bug corrige)
    F = (2 * rand + 1) * (1 - t / Tmax);

    for i = 1:Npop
        % Selection du vautour de reference selon p1
        if rand < p1
            R = Best1_pos;
        else
            R = Best2_pos;
        end

        if abs(F) >= 1
            % ============================================================
            % EXPLORATION — deux strategies selon p2
            % ============================================================
            if rand < p2
                % Strategie 1 : deplacement vers reference avec perturbation
                A = 2 * rand(1,dim) - 1;
                X(i,:) = R - abs(A .* R - X(i,:)) .* F;
            else
                % Strategie 2 : deplacement aleatoire dans l'espace
                X(i,:) = R - F + rand(1,dim) .* ((ub - lb) + lb);
            end

        else
            % ============================================================
            % EXPLOITATION — deux phases selon |F|
            % ============================================================
            if abs(F) >= 0.5
                % Phase 1 : rotation agressive (Rotational Flight)
                if rand < p3
                    % Moyenne ponderee des deux meilleurs
                    D1 = abs(X(i,:) - Best1_pos);
                    D2 = abs(X(i,:) - Best2_pos);
                    X(i,:) = (alpha * Best1_pos + (1-alpha) * Best2_pos) ...
                              - F * (D1 + D2) / 2;
                else
                    X(i,:) = R - abs(F) * rand(1,dim) .* (R - X(i,:));
                end
            else
                % Phase 2 : exploitation douce (Siege Tendency)
                if rand < p3
                    % Levy flight vers Best1
                    LF = levy_flight(1, dim, 1.5);
                    S1 = Best1_pos - (Best1_pos .* X(i,:)) / ...
                         (Best1_pos - X(i,:) + 1e-10);
                    S2 = Best2_pos - (Best2_pos .* X(i,:)) / ...
                         (Best2_pos - X(i,:) + 1e-10);
                    X(i,:) = Best1_pos - (S1 + S2) .* LF;
                else
                    % Mouvement spiral
                    theta = rand * 2 * pi;
                    r1    = rand; r2 = rand;
                    X(i,:) = abs(R - X(i,:)) .* cos(2*pi*theta) .* r1 + ...
                              R + sin(2*pi*theta) .* r2 .* (Best1_pos - X(i,:));
                end
            end
        end

        % Clampage aux bornes
        X(i,:) = max(lb, min(ub, X(i,:)));

        % Evaluation
        f = fobj(X(i,:));
        fit(i) = f;

        % Mise a jour hierarchie des deux meilleurs
        if f < Best1_fit
            Best2_fit = Best1_fit;  Best2_pos = Best1_pos;
            Best1_fit = f;          Best1_pos = X(i,:);
        elseif f < Best2_fit
            Best2_fit = f;  Best2_pos = X(i,:);
        end
    end

    Convergence_curve(t) = Best1_fit;
end

Best_fit = Best1_fit;
Best_pos = Best1_pos;

end
