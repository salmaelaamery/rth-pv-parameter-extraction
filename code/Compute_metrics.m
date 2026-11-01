function [RMSE, MAE, MBE, SSE] = Compute_metrics(I_exp, I_model)
% =========================================================================
% Compute_metrics.m
% =========================================================================
% Calcule les metriques d'erreur entre les courants experimentaux et estimes.
%
% Entrees :
%   I_exp   : vecteur des courants experimentaux [A]
%   I_model : vecteur des courants estimes par le modele [A]
%
% Sorties :
%   RMSE : Root Mean Square Error (racine de l'erreur quadratique moyenne)
%          RMSE = sqrt( (1/N) * sum((I_exp - I_model)^2) )
%
%   MAE  : Mean Absolute Error (erreur absolue moyenne)
%          MAE  = (1/N) * sum(|I_exp - I_model|)
%
%   MBE  : Mean Bias Error (erreur de biais moyenne)
%          MBE  = (1/N) * sum(I_model - I_exp)
%          > 0 : surestimation | < 0 : sous-estimation
%
%   SSE  : Sum of Squared Errors (somme des carres des erreurs)
%          SSE  = sum((I_exp - I_model)^2)
%
% Auteur : adapte pour l'extraction des parametres PV - RTH
% =========================================================================

% --- S'assurer que les vecteurs sont en colonnes -------------------------
I_exp   = I_exp(:);
I_model = I_model(:);

N = length(I_exp);

% --- Calcul des erreurs --------------------------------------------------
residuals = I_exp - I_model;    % erreurs signees : exp - modele

% --- RMSE ----------------------------------------------------------------
RMSE = sqrt(sum(residuals.^2) / N);

% --- MAE -----------------------------------------------------------------
MAE  = sum(abs(residuals)) / N;

% --- MBE -----------------------------------------------------------------
% Convention : MBE positif signifie que le modele sous-estime
MBE  = sum(I_model - I_exp) / N;

% --- SSE -----------------------------------------------------------------
SSE  = sum(residuals.^2);

end
