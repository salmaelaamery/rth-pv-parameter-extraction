function I_model = PV_current_model(x, V_exp, Ns, T)
% =========================================================================
% PV_current_model.m
% =========================================================================
% Calcule le courant du modele single-diode pour chaque tension V_exp.
%
% Modele implicite :
%   I = Iph - I0*(exp((V + I*Rs)/(n*Vt)) - 1) - (V + I*Rs)/Rsh
%
% Resolution par la methode de Newton-Raphson vectorisee (robuste).
%
% Entrees :
%   x     : vecteur de parametres [Iph, I0, Rs, Rsh, n]
%   V_exp : vecteur des tensions experimentales (colonne, en V)
%   Ns    : nombre de cellules en serie
%   T     : temperature en Kelvin
%
% Sortie :
%   I_model : vecteur des courants estimes (meme taille que V_exp)
%
% Auteur : adapte pour l'extraction des parametres PV - RTH
% =========================================================================

% --- Constantes physiques -------------------------------------------------
k = 1.380649e-23;     % Constante de Boltzmann [J/K]
q = 1.602176634e-19;  % Charge elementaire [C]

% --- Extraction des parametres -------------------------------------------
Iph = x(1);   % Courant photo-genere [A]
I0  = x(2);   % Courant de saturation inverse [A]
Rs  = x(3);   % Resistance serie [Ohm]
Rsh = x(4);   % Resistance shunt [Ohm]
n   = x(5);   % Facteur d'idealite [-]

% --- Tension thermique ---------------------------------------------------
Vt = Ns * k * T / q;   % [V]

% --- Protection contre les valeurs numeriquement dangereuses -------------
% Eviter division par zero
if Rsh < 1e-12
    Rsh = 1e-12;
end
if abs(n * Vt) < 1e-12
    n = 1e-10 / Vt;
end

% --- Initialisation de la solution (approximation : I ≈ Iph) ------------
V_exp   = V_exp(:);      % S'assurer que c'est un vecteur colonne
N_pts   = length(V_exp);
I_model = Iph * ones(N_pts, 1);

% --- Parametres de Newton-Raphson ----------------------------------------
max_iter  = 200;    % Nombre maximum d'iterations
tol       = 1e-10;  % Tolerance de convergence
exp_limit = 700;    % Limite pour eviter overflow de l'exponentielle

% --- Boucle Newton-Raphson -----------------------------------------------
for iter = 1:max_iter

    % Argument de l'exponentielle (clippe pour la stabilite numerique)
    arg = (V_exp + I_model .* Rs) ./ (n * Vt);
    arg = min(arg, exp_limit);   % Eviter exp(+Inf)

    exp_val = exp(arg);

    % f(I) = I - Iph + I0*(exp_val - 1) + (V + I*Rs)/Rsh
    f  = I_model - Iph ...
         + I0 .* (exp_val - 1) ...
         + (V_exp + I_model .* Rs) ./ Rsh;

    % f'(I) = 1 + I0*Rs/(n*Vt)*exp_val + Rs/Rsh
    df = 1 ...
         + I0 .* Rs ./ (n * Vt) .* exp_val ...
         + Rs ./ Rsh;

    % Protection contre division par zero dans la derivee
    df(abs(df) < 1e-15) = 1e-15;

    % Mise a jour Newton-Raphson
    delta = f ./ df;
    I_model = I_model - delta;

    % Convergence atteinte ?
    if max(abs(delta)) < tol
        break;
    end
end

% --- Remplacer les valeurs non-physiques (NaN, Inf, negatifs extremes) ---
I_model(~isfinite(I_model)) = 0;

end
