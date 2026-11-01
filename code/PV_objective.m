function fitness = PV_objective(x, V_exp, I_exp, Ns, T)
% =========================================================================
% PV_objective.m
% =========================================================================
% Fonction objectif pour l'extraction des parametres du modele PV
% single-diode par algorithme metaheuristique (RTH, PSO, GWO, etc.).
%
% Minimise le RMSE entre courants experimentaux et courants estimes.
%
% CONVENTION DE PARAMETRES (espace d'optimisation) :
%   x(1) = Iph           [A]         courant photo-genere
%   x(2) = log10(I0)     [log10(A)]  LOG du courant de saturation
%   x(3) = Rs            [Ohm]       resistance serie
%   x(4) = Rsh           [Ohm]       resistance shunt
%   x(5) = n             [-]         facteur d'idealite
%
%   La transformation log10 sur I0 est INDISPENSABLE : I0 varie sur
%   9 ordres de grandeur (1e-12 a 1e-3 A). Sans cette transformation,
%   l'algorithme explore l'espace en echelle lineaire et ne trouve jamais
%   la region autour de ~1e-7 A caracteristique d'une cellule silicium.
%
% Utilisation avec RTH :
%   fobj = @(x) PV_objective(x, V_exp, I_exp, Ns, T);
%   [Cost, Pos, Conv] = RTH(Npop, Tmax, lb, ub, dim, fobj);
%   Bornes : lb = [0, -12, 0,   1,   1]
%            ub = [2,  -3, 1, 2000,  2]
%
% Auteur : adapte pour l'extraction des parametres PV - RTH
% =========================================================================

% --- Conversion espace optimisation -> parametres physiques --------------
% x(2) = log10(I0)  =>  I0 = 10^x(2)
x_phys    = x;
x_phys(2) = 10^(x(2));   % Retour en echelle lineaire

% --- Verification validite physique minimale -----------------------------
if x_phys(4) <= 0 || x_phys(5) <= 0   % Rsh > 0 et n > 0
    fitness = 1e10;
    return;
end
if x_phys(1) < 0 || x_phys(3) < 0    % Iph >= 0 et Rs >= 0
    fitness = 1e10;
    return;
end

% --- Calcul du courant estime (avec parametres physiques) ----------------
try
    I_model = PV_current_model(x_phys, V_exp, Ns, T);
catch
    fitness = 1e10;
    return;
end

% --- Verification de la validite numerique --------------------------------
if any(~isfinite(I_model))
    fitness = 1e10;
    return;
end

% --- Calcul du RMSE -------------------------------------------------------
N = length(I_exp);
residuals = I_exp(:) - I_model(:);
fitness   = sqrt(sum(residuals.^2) / N);

% Penalite si la valeur est anormale
if ~isfinite(fitness)
    fitness = 1e10;
end

end
