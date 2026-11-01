function [Best_fit, Best_pos, Convergence_curve] = HHO(Npop, Tmax, lb, ub, dim, fobj)
% =========================================================================
% HHO.m  —  Harris Hawks Optimization
% =========================================================================
% Heidari et al. (2019), Future Generation Computer Systems.
%
% Interface standard :
%   [Best_fit, Best_pos, Convergence_curve] = HHO(Npop, Tmax, lb, ub, dim, fobj)
%
% Phases :
%   - Exploration : sauts aleatoires vers lapin ou position moyenne
%   - Transition  : E (energie d'echappement) decroit avec le temps
%   - Exploitation : attaque douce / dure, avec/sans rafale d'energie
%   Levy flight pour la phase d'attaque
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

Convergence_curve = zeros(1, Tmax);

% --- Fonction Levy -------------------------------------------------------
    function L = levy(n, m, beta)
        num   = gamma(1 + beta) * sin(pi * beta / 2);
        den   = gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2);
        sigma = (num / den)^(1 / beta);
        u     = randn(n, m) * sigma;
        v     = randn(n, m);
        L     = u ./ (abs(v).^(1 / beta));
    end

% --- Boucle principale ---------------------------------------------------
for t = 1:Tmax
    % Energie d'echappement du lapin
    E1 = 2 * (1 - t / Tmax);   % decroit de 2 a 0

    % Mise a jour position moyenne
    X_mean = mean(X, 1);

    for i = 1:Npop
        E0  = 2 * rand - 1;      % energie initiale dans [-1, 1]
        Esc = E1 * E0;            % energie d'echappement

        if abs(Esc) >= 1
            % === EXPLORATION ===
            q = rand;
            rand_idx = randi(Npop);
            X_rand = X(rand_idx, :);

            if q >= 0.5
                % Saut vers un individu aleatoire
                X(i,:) = X_rand - rand(1,dim) .* abs(X_rand - 2*rand(1,dim).*Best_pos);
            else
                % Saut vers lapin via position moyenne
                X(i,:) = (Best_pos - X_mean) - rand(1,dim) .* ...
                          (lb + rand(1,dim) .* (ub - lb));
            end

        else
            % === EXPLOITATION ===
            r = rand;
            delta = Best_pos - X(i,:);

            if r >= 0.5 && abs(Esc) >= 0.5
                % Attaque douce
                X(i,:) = delta - Esc .* abs(delta);

            elseif r >= 0.5 && abs(Esc) < 0.5
                % Attaque dure
                X(i,:) = Best_pos - Esc .* abs(delta);

            elseif r < 0.5 && abs(Esc) >= 0.5
                % Attaque douce avec Levy
                LF = levy(1, dim, 1.5);
                Y  = Best_pos - Esc .* abs(delta);
                Y  = max(lb, min(ub, Y));
                Z  = Y + rand(1,dim) .* LF;
                Z  = max(lb, min(ub, Z));

                fY = fobj(Y);
                fZ = fobj(Z);
                if fY < fit(i)
                    X(i,:) = Y; fit(i) = fY;
                end
                if fZ < fit(i)
                    X(i,:) = Z; fit(i) = fZ;
                end
                % skip normal update below for this branch
                if fY < Best_fit
                    Best_fit = fY; Best_pos = Y;
                end
                if fZ < Best_fit
                    Best_fit = fZ; Best_pos = Z;
                end
                continue;  % passe a l'individu suivant

            else
                % Attaque dure avec Levy
                LF = levy(1, dim, 1.5);
                Y  = Best_pos - Esc .* abs(delta);
                Y  = max(lb, min(ub, Y));
                Z  = Y + rand(1,dim) .* LF;
                Z  = max(lb, min(ub, Z));

                fY = fobj(Y);
                fZ = fobj(Z);
                if fY < fit(i)
                    X(i,:) = Y; fit(i) = fY;
                end
                if fZ < fit(i)
                    X(i,:) = Z; fit(i) = fZ;
                end
                if fY < Best_fit
                    Best_fit = fY; Best_pos = Y;
                end
                if fZ < Best_fit
                    Best_fit = fZ; Best_pos = Z;
                end
                continue;
            end
        end

        % Clampage & evaluation (branches sans Levy)
        X(i,:) = max(lb, min(ub, X(i,:)));
        f = fobj(X(i,:));
        fit(i) = f;

        if f < Best_fit
            Best_fit = f;
            Best_pos = X(i,:);
        end
    end

    Convergence_curve(t) = Best_fit;
end

end
