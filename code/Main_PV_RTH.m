% =========================================================================
% Main_PV_RTH.m
% =========================================================================
% Extraction des parametres du modele PV single-diode
% par le Red-Tailed Hawk Algorithm (RTH)
%
% Modele : I = Iph - I0*[exp((V+I*Rs)/(n*Vt))-1] - (V+I*Rs)/Rsh
% Parametres : x = [Iph, I0, Rs, Rsh, n]
%
% Reference : RTH - DOI: 10.1038/s41598-023-38778-3
%
% Datasets benchmark supportes :
%   - Shell SM55 55W monocristallin (defaut - dataset Chaibi 2018/2019)
%   - RTC France 57mm (benchmark classique litterature)
%   - Photowatt PWP201
%   - Chargement depuis fichier CSV (deux colonnes : V, I)
%
% Pour lancer : executer ce script dans MATLAB.
%   Les fichiers requis dans le meme repertoire :
%     RTH.m | PV_objective.m | PV_current_model.m | Compute_metrics.m
%
% Auteur : adapte pour l'extraction des parametres PV - RTH
% =========================================================================

clear all; clc; close all;

%% ========================================================================
%% SECTION 1 : CHOIX DU DATASET
%% ========================================================================
% Selectionner le dataset en decommentant la ligne souhaitee :

% DATASET = 'RTC_France';         % Cellule Si 57mm - benchmark classique
DATASET = 'SM55';                  % <-- Shell SM55 55W monocristallin (ACTIF)
% DATASET = 'Photowatt_PWP201';   % Module polycristallin 36 cellules
% DATASET = 'CSV';                 % Charger depuis fichier CSV

switch DATASET

    % ------------------------------------------------------------------
    case 'RTC_France'
    % RTC France - Cellule silicium 57mm
    % Mesures a T = 33 degC, Ns = 1 cellule
    % Source : Easwarakhanthan et al. (1986)
    % ------------------------------------------------------------------
        T  = 33 + 273.15;   % Temperature [K]
        Ns = 1;             % Cellules en serie

        V_exp = [-0.2057; -0.1291; -0.0588;  0.0057;  0.0646;  0.1185;
                  0.1678;  0.2132;  0.2545;  0.2924;  0.3269;  0.3585;
                  0.3873;  0.4137;  0.4373;  0.4590;  0.4784;  0.4960;
                  0.5119;  0.5265;  0.5398;  0.5521;  0.5633;  0.5736;
                  0.5833;  0.5900];

        I_exp = [ 0.7640;  0.7620;  0.7605;  0.7605;  0.7600;  0.7590;
                  0.7570;  0.7570;  0.7555;  0.7540;  0.7505;  0.7465;
                  0.7385;  0.7280;  0.7065;  0.6755;  0.6320;  0.5730;
                  0.4990;  0.4130;  0.3165;  0.2120;  0.1035; -0.0100;
                 -0.1230; -0.2100];

    % ------------------------------------------------------------------
    case 'Photowatt_PWP201'
    % Photowatt PWP201 - Module polysilicium, 36 cellules
    % Mesures a T = 45 degC, Irradiance = 1000 W/m2
    % Source : Sera et al. (2007)
    % ------------------------------------------------------------------
        T  = 45 + 273.15;
        Ns = 36;

        V_exp = [0.1248; 1.8093; 3.3511; 4.7622; 6.0538; 7.2364;
                 8.3189; 9.3097; 10.2163; 11.0449; 11.8018; 12.4929;
                 13.1231; 13.6983; 14.2232; 14.6995; 15.1328; 15.5240;
                 15.8734; 16.1795; 16.4416; 16.6576];

        I_exp = [1.0315; 1.0300; 1.0260; 1.0200; 1.0100; 0.9960;
                 0.9760; 0.9510; 0.9200; 0.8830; 0.8390; 0.7870;
                 0.7270; 0.6580; 0.5800; 0.4910; 0.3930; 0.2870;
                 0.1740; 0.0540; -0.0630; -0.1300];

    % ------------------------------------------------------------------
    case 'SM55'
    % Shell SM55 - Module monocristallin silicium, 36 cellules
    % Conditions STC : T = 25 degC, lambda = 1000 W/m2
    % Source : Chaibi et al. (Solar Energy 2018 + Data in Brief 2019)
    %          Courbe I-V modellee avec les parametres extraits par Chaibi :
    %          Rs=0.1124 Ohm, Ior=4.8424e-6 A, gamma=1.7411, Rsh=6500 Ohm
    % Specifications constructeur (STC) :
    %   Pm=55W, Vm=17.4V, Im=3.15A, Voc=21.7V, Isc=3.45A, Ns=36
    %
    % Valeurs de reference Chaibi (2018) pour validation :
    %   Rs    = 0.1124 Ohm
    %   I0    = 4.8424e-6 A  (Ior)
    %   n     = 1.7411       (gamma)
    %   Rsh   = 6500 Ohm     (mesure en laboratoire)
    %   Iph   ~ 3.450 A
    % ------------------------------------------------------------------
        T  = 25 + 273.15;   % Temperature [K] - STC
        Ns = 36;            % Cellules en serie

        % Points V-I modelises depuis les parametres Chaibi (2018)
        % 25 points bien distribues de 0 V a la zone proche de Voc (~21.65 V)
        % Distribution dense dans la zone du coude (19-21.65 V)
        V_exp = [ 0.0000;  1.3571;  2.7143;  4.0714;  5.4286;
                  6.7857;  8.1429;  9.5000; 10.8571; 12.2143;
                 13.5714; 14.9286; 16.2857; 17.6429; 19.0000;
                 19.2650; 19.5300; 19.7950; 20.0600; 20.3250;
                 20.5900; 20.8550; 21.1200; 21.3850; 21.6500];

        I_exp = [ 3.4500;  3.4498;  3.4496;  3.4493;  3.4490;
                  3.4485;  3.4478;  3.4463;  3.4431;  3.4360;
                  3.4198;  3.3826;  3.2972;  3.1029;  2.6709;
                  2.5403;  2.3893;  2.2152;  2.0150;  1.7856;
                  1.5239;  1.2264;  0.8900;  0.5112;  0.0871];

    % ------------------------------------------------------------------
    case 'CSV'
    % Chargement depuis un fichier CSV (deux colonnes : V et I)
    % Modifier le nom du fichier ci-dessous
    % ------------------------------------------------------------------
        T  = 25 + 273.15;   % <-- Modifier la temperature
        Ns = 1;             % <-- Modifier le nombre de cellules

        csv_file = 'data_PV.csv';   % <-- Nom du fichier CSV
        if ~isfile(csv_file)
            error('Fichier CSV introuvable : %s', csv_file);
        end
        data  = readmatrix(csv_file);
        V_exp = data(:, 1);
        I_exp = data(:, 2);
        fprintf('Donnees chargees depuis : %s (%d points)\n', csv_file, length(V_exp));

    otherwise
        error('Dataset inconnu : %s', DATASET);
end

% S'assurer que les vecteurs sont en colonnes
V_exp = V_exp(:);
I_exp = I_exp(:);

if strcmp(DATASET, 'SM55')
    data_label = 'Reference points';
else
    data_label = 'Experimental data';
end

%% ========================================================================
%% SECTION 2 : CONSTANTES PHYSIQUES (reference uniquement)
%% ========================================================================
k = 1.380649e-23;     % Constante de Boltzmann [J/K]
q = 1.602176634e-19;  % Charge elementaire [C]
% Vt = Ns * k * T / q  --> calcule dans PV_current_model.m

%% ========================================================================
%% SECTION 3 : PARAMETRES D'OPTIMISATION
%% ========================================================================
dim = 5;    % Dimension : [Iph, I0, Rs, Rsh, n]

% --- Bornes des parametres (espace d'optimisation) -----------------------
%
%  IMPORTANT : x(2) = log10(I0), pas I0 directement.
%  Cette transformation est cle pour la convergence : I0 varie sur
%  9 ordres de grandeur. En echelle log, l'algorithme explore
%  uniformement toute la plage [1e-12, 1e-3] A.
%
%                  Iph[A]  log10(I0)  Rs[Ohm]  Rsh[Ohm]  n[-]
% SM55             [0, 5]  [-12, -3]  [0, 0.5] [1, 8000] [1, 2]
% RTC France       [0, 1]  [-12, -3]  [0, 0.5] [1, 100]  [1, 2]
% Photowatt PWP201 [0, 2]  [-12, -3]  [0, 0.5] [1, 2000] [1, 2]

switch DATASET
    case 'SM55'
        %       Iph[A]  log10(I0)   Rs[Ohm]  Rsh[Ohm]  n[-]
        lb = [  0,      -12,        0,        1,         1  ];
        ub = [  5,      -3,         0.5,      8000,      2  ];
    case 'RTC_France'
        lb = [  0,      -12,        0,        1,         1  ];
        ub = [  1,      -3,         0.5,      100,       2  ];
    case 'Photowatt_PWP201'
        lb = [  0,      -12,        0,        1,         1  ];
        ub = [  2,      -3,         0.5,      2000,      2  ];
    otherwise
        lb = [  0,      -12,        0,        1,         1  ];
        ub = [  10,     -3,         1,        5000,      2  ];
end
%
% Les parametres physiques sont recuperes par : I0 = 10^(x(2))

% --- Parametres RTH ------------------------------------------------------
Npop  = 50;    % Taille de la population (augmentee pour meilleure exploration)
Tmax  = 1000;  % Nombre maximum d'iterations
Nruns = 30;    % Nombre de runs independants

% --- Fonction objectif ---------------------------------------------------
fobj = @(x) PV_objective(x, V_exp, I_exp, Ns, T);

%% ========================================================================
%% SECTION 4 : EXECUTIONS INDEPENDANTES (30 RUNS)
%% ========================================================================
% Pre-allocation des structures de resultats
results_fitness = zeros(Nruns, 1);
results_params  = zeros(Nruns, dim);
results_RMSE    = zeros(Nruns, 1);
results_MAE     = zeros(Nruns, 1);
results_MBE     = zeros(Nruns, 1);
results_SSE     = zeros(Nruns, 1);
results_time    = zeros(Nruns, 1);
conv_curves     = zeros(Nruns, Tmax);

fprintf('\n');
fprintf('=================================================================\n');
fprintf('   PV Parameter Extraction - Red-Tailed Hawk Algorithm (RTH)\n');
fprintf('=================================================================\n');
fprintf('   Dataset  : %s\n', DATASET);
fprintf('   T        : %.2f K (%.2f degC)\n', T, T - 273.15);
fprintf('   Ns       : %d cellule(s)\n', Ns);
fprintf('   N points : %d\n', length(V_exp));
fprintf('   Dim      : %d   |  Npop : %d   |  Tmax : %d\n', dim, Npop, Tmax);
fprintf('   Nruns    : %d\n', Nruns);
fprintf('-----------------------------------------------------------------\n');
fprintf('   Run  |  RMSE            |  Temps (s)\n');
fprintf('-----------------------------------------------------------------\n');

rng('shuffle');   % Graine aleatoire differente a chaque execution

for run = 1:Nruns
    tic;
    [Best_fit, Best_pos, Convergence_curve] = RTH(Npop, Tmax, lb, ub, dim, fobj);
    results_time(run) = toc;

    % Stocker les resultats bruts
    results_fitness(run)   = Best_fit;
    results_params(run, :) = Best_pos;
    conv_curves(run, :)    = Convergence_curve;

    % Calculer les metriques pour ce run
    % Conversion log10(I0) -> I0 physique avant appel au modele
    Best_pos_phys    = Best_pos;
    Best_pos_phys(2) = 10^(Best_pos(2));
    I_model = PV_current_model(Best_pos_phys, V_exp, Ns, T);
    [results_RMSE(run), results_MAE(run), results_MBE(run), results_SSE(run)] = ...
        Compute_metrics(I_exp, I_model);

    fprintf('   %3d  |  %.6e      |  %.2f s\n', run, results_RMSE(run), results_time(run));
end

%% ========================================================================
%% SECTION 5 : IDENTIFICATION DU MEILLEUR RUN
%% ========================================================================
[~, best_run_idx] = min(results_RMSE);

% Parametres optimises (espace log pour I0)
best_params = results_params(best_run_idx, :);

% Parametres physiques (I0 en lineaire) pour affichage et tracé
best_params_phys    = best_params;
best_params_phys(2) = 10^(best_params(2));

% Parametres physiques de tous les runs (pour affichage/export)
display_params = results_params;
display_params(:,2) = 10.^(results_params(:,2));

best_I_model = PV_current_model(best_params_phys, V_exp, Ns, T);
[RMSE_best, MAE_best, MBE_best, SSE_best] = Compute_metrics(I_exp, best_I_model);

%% ========================================================================
%% SECTION 6 : TABLEAU DETAILLE DES 30 RUNS
%% ========================================================================
fprintf('\n');
fprintf('=================================================================\n');
fprintf('   TABLEAU DETAILLE DES %d RUNS\n', Nruns);
fprintf('=================================================================\n');

header = sprintf('%-4s | %-12s | %-9s | %-11s | %-9s | %-9s | %-7s | %-12s | %-12s | %-12s | %-12s | %-8s', ...
    'Run', 'Fitness', 'Iph (A)', 'I0 (A)', 'Rs (Ohm)', 'Rsh (Ohm)', 'n', ...
    'RMSE', 'MAE', 'MBE', 'SSE', 'Temps(s)');
fprintf('%s\n', header);
fprintf('%s\n', repmat('-', 1, length(header)));

for run = 1:Nruns
    marker = '';
    if run == best_run_idx
        marker = ' <-- BEST';
    end
    % display_params(:,2) contient deja I0 en lineaire (10^x(2))
    fprintf('%-4d | %-12.4e | %-9.6f | %-11.4e | %-9.6f | %-9.4f | %-7.5f | %-12.4e | %-12.4e | %-12.4e | %-12.4e | %-8.2f%s\n', ...
        run, results_fitness(run), ...
        display_params(run,1), display_params(run,2), ...
        display_params(run,3), display_params(run,4), display_params(run,5), ...
        results_RMSE(run), results_MAE(run), results_MBE(run), results_SSE(run), ...
        results_time(run), marker);
end

%% ========================================================================
%% SECTION 7 : RESUME STATISTIQUE
%% ========================================================================
fprintf('\n');
fprintf('=================================================================\n');
fprintf('   RESUME STATISTIQUE SUR %d RUNS\n', Nruns);
fprintf('=================================================================\n');
fprintf('   Best RMSE  : %.6e\n', min(results_RMSE));
fprintf('   Mean RMSE  : %.6e\n', mean(results_RMSE));
fprintf('   Std RMSE   : %.6e\n', std(results_RMSE));
fprintf('   Worst RMSE : %.6e\n', max(results_RMSE));
fprintf('-----------------------------------------------------------------\n');
fprintf('   Meilleurs parametres (Run %d) :\n', best_run_idx);
fprintf('     Iph = %.8f  A\n',   best_params_phys(1));
fprintf('     I0  = %.8e  A\n',   best_params_phys(2));
fprintf('     Rs  = %.8f  Ohm\n', best_params_phys(3));
fprintf('     Rsh = %.6f  Ohm\n', best_params_phys(4));
fprintf('     n   = %.8f\n',      best_params_phys(5));
fprintf('-----------------------------------------------------------------\n');
fprintf('   Metriques (meilleur run) :\n');
fprintf('     RMSE = %.6e\n', RMSE_best);
fprintf('     MAE  = %.6e\n', MAE_best);
fprintf('     MBE  = %.6e\n', MBE_best);
fprintf('     SSE  = %.6e\n', SSE_best);
fprintf('=================================================================\n');

%% ========================================================================
%% SECTION 8 : FIGURES
%% ========================================================================
% Parametres graphiques globaux (style article scientifique)
set(0, 'DefaultAxesFontName',  'Times New Roman');
set(0, 'DefaultAxesFontSize',  13);
set(0, 'DefaultTextFontName',  'Times New Roman');
set(0, 'DefaultLineLineWidth', 1.8);

color_exp   = [0.10, 0.10, 0.10];   % noir pour les points de reference/mesure
color_model = [0.85, 0.10, 0.10];   % rouge pour le modele RTH
color_conv  = [0.10, 0.40, 0.75];   % bleu pour la convergence

% ---------------------------------------------------------
% Figure 1 : Courbe I-V (experimental vs modele RTH)
% ---------------------------------------------------------
fig1 = figure('Position', [80, 80, 680, 480], 'Color', 'w');
hold on;
plot(V_exp, I_exp, 'o', ...
    'Color', color_exp, 'MarkerFaceColor', color_exp, ...
    'MarkerSize', 6, 'DisplayName', data_label);
plot(V_exp, best_I_model, '-', ...
    'Color', color_model, 'LineWidth', 2, 'DisplayName', 'RTH model');
hold off;
xlabel('Voltage  V  (V)', 'FontWeight', 'bold');
ylabel('Current  I  (A)', 'FontWeight', 'bold');
title(['I' char(8211) 'V Characteristic: ' strrep(DATASET,'_',' ')], ...
      'FontWeight', 'bold', 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 12, 'Box', 'on');
grid on; box on;
set(gca, 'GridAlpha', 0.25, 'GridLineStyle', '--');

% ---------------------------------------------------------
% Figure 2 : Courbe P-V (experimental vs modele RTH)
% ---------------------------------------------------------
P_exp   = V_exp .* I_exp;
P_model = V_exp .* best_I_model;

fig2 = figure('Position', [140, 140, 680, 480], 'Color', 'w');
hold on;
plot(V_exp, P_exp, 'o', ...
    'Color', color_exp, 'MarkerFaceColor', color_exp, ...
    'MarkerSize', 6, 'DisplayName', data_label);
plot(V_exp, P_model, '-', ...
    'Color', color_model, 'LineWidth', 2, 'DisplayName', 'RTH model');
hold off;
xlabel('Voltage  V  (V)', 'FontWeight', 'bold');
ylabel('Power  P  (W)', 'FontWeight', 'bold');
title(['P' char(8211) 'V Characteristic: ' strrep(DATASET,'_',' ')], ...
      'FontWeight', 'bold', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 12, 'Box', 'on');
grid on; box on;
set(gca, 'GridAlpha', 0.25, 'GridLineStyle', '--');

% ---------------------------------------------------------
% Figure 3 : Courbe de convergence (meilleur run)
% ---------------------------------------------------------
fig3 = figure('Position', [200, 200, 680, 430], 'Color', 'w');
semilogy(1:Tmax, conv_curves(best_run_idx, :), '-', ...
    'Color', color_conv, 'LineWidth', 2);
xlabel('Iteration', 'FontWeight', 'bold');
ylabel('Best RMSE (log scale)', 'FontWeight', 'bold');
title(sprintf('Convergence Curve \x2013 Best Run (Run %d)', best_run_idx), ...
      'FontWeight', 'bold', 'FontSize', 14);
legend('RTH', 'Location', 'northeast', 'FontSize', 12, 'Box', 'on');
grid on; box on;
set(gca, 'GridAlpha', 0.25, 'GridLineStyle', '--');

% ---------------------------------------------------------
% Figure 4 : Boxplot des RMSE sur 30 runs
% ---------------------------------------------------------
fig4 = figure('Position', [260, 260, 520, 480], 'Color', 'w');
boxplot(results_RMSE, 'Labels', {'RTH'}, ...
        'Widths', 0.45, 'Colors', color_model, ...
        'Symbol', 'r+', 'OutlierSize', 7);
ylabel('RMSE', 'FontWeight', 'bold');
title(sprintf('RMSE Distribution over %d Runs', Nruns), ...
      'FontWeight', 'bold', 'FontSize', 14);
grid on; box on;
set(gca, 'GridAlpha', 0.25, 'GridLineStyle', '--');

% Annoter le graphique avec les statistiques
hold on;
yline(mean(results_RMSE), '--b', sprintf('Mean=%.3e', mean(results_RMSE)), ...
      'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.2);
hold off;

%% ========================================================================
%% SECTION 9 : EXPORT DES RESULTATS VERS EXCEL
%% ========================================================================
filename_xlsx = sprintf('RTH_PV_results_%s.xlsx', DATASET);

% --- Feuille 1 : Tous les runs -------------------------------------------
run_col   = (1:Nruns)';
% display_params(:,2) = I0 en lineaire (A), les autres en valeurs directes
T_runs    = table(run_col, results_fitness, ...
    display_params(:,1), display_params(:,2), ...
    display_params(:,3), display_params(:,4), display_params(:,5), ...
    results_RMSE, results_MAE, results_MBE, results_SSE, results_time, ...
    'VariableNames', {'Run', 'Best_Fitness', ...
    'Iph_A', 'I0_A', 'Rs_Ohm', 'Rsh_Ohm', 'n', ...
    'RMSE', 'MAE', 'MBE', 'SSE', 'Time_s'});

writetable(T_runs, filename_xlsx, 'Sheet', 'All_Runs');

% --- Feuille 2 : Meilleurs parametres ------------------------------------
T_best = table(best_run_idx, best_params_phys(1), best_params_phys(2), ...
    best_params_phys(3), best_params_phys(4), best_params_phys(5), ...
    RMSE_best, MAE_best, MBE_best, SSE_best, ...
    'VariableNames', {'Best_Run', 'Iph_A', 'I0_A', 'Rs_Ohm', 'Rsh_Ohm', 'n', ...
    'RMSE', 'MAE', 'MBE', 'SSE'});

writetable(T_best, filename_xlsx, 'Sheet', 'Best_Parameters');

% --- Feuille 3 : Statistiques --------------------------------------------
T_stats = table(min(results_RMSE), mean(results_RMSE), ...
    std(results_RMSE), max(results_RMSE), ...
    sum(results_time), ...
    'VariableNames', {'Best_RMSE', 'Mean_RMSE', 'Std_RMSE', 'Worst_RMSE', ...
    'Total_Time_s'});

writetable(T_stats, filename_xlsx, 'Sheet', 'Statistics');

% --- Feuille 4 : Donnees experimentales vs modele ------------------------
T_IV = table(V_exp, I_exp, best_I_model, ...
    I_exp - best_I_model, ...
    V_exp .* I_exp, V_exp .* best_I_model, ...
    'VariableNames', {'V_exp_V', 'I_exp_A', 'I_model_A', ...
    'Error_A', 'P_exp_W', 'P_model_W'});

writetable(T_IV, filename_xlsx, 'Sheet', 'IV_Comparison');

fprintf('\n   Resultats exportes vers : %s\n', filename_xlsx);

% --- Sauvegarder les figures en haute resolution -------------------------
saveas(fig1, sprintf('IV_curve_%s.png', DATASET));
saveas(fig2, sprintf('PV_curve_%s.png', DATASET));
saveas(fig3, sprintf('Convergence_%s.png', DATASET));
saveas(fig4, sprintf('RMSE_boxplot_%s.png', DATASET));

% Pour article : export en PDF ou EMF (vecteur)
% print(fig1, sprintf('IV_curve_%s', DATASET), '-dpdf', '-r300');

fprintf('   Figures sauvegardees (PNG).\n');
fprintf('=================================================================\n');
fprintf('   Termine avec succes.\n');
fprintf('=================================================================\n\n');

%% ========================================================================
%% NOTE : Pour comparer plusieurs algorithmes
%% ========================================================================
% Il suffit de remplacer la ligne :
%   [Best_fit, Best_pos, Convergence_curve] = RTH(Npop, Tmax, lb, ub, dim, fobj);
% par l'appel d'un autre algorithme (PSO, GWO, WOA, HHO, AVOA, MGO...)
% avec la meme signature : [Best_fit, Best_pos, Convergence_curve] = ALGO(...)
% Toutes les fonctions objectif et de metriques restent identiques.
