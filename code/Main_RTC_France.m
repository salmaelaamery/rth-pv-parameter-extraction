% =========================================================================
% Main_RTC_France.m  —  RTC France 57mm Solar Cell Benchmark
% =========================================================================
% Multi-algorithm comparison for single-diode PV parameter extraction.
% Cell   : RTC France 57 mm silicon solar cell (Ns=1, T=33 degC)
% Dataset: Easwarakhanthan et al., Int. J. Solar Energy, 4(1):1-12, 1986
%
% Algorithms compared (identical conditions to SM55 benchmark):
%   RTH  — Red-Tailed Hawk Algorithm (Ferahtia et al. 2023)
%   PSO  — Particle Swarm Optimization
%   GWO  — Grey Wolf Optimizer
%   WOA  — Whale Optimization Algorithm
%   HHO  — Harris Hawks Optimization
%   AVOA — African Vultures Optimization Algorithm
%   MGO  — Mountain Gazelle Optimizer
%
% Settings (same as SM55 benchmark for fair cross-dataset comparison):
%   Npop = 50,  Tmax = 1000,  Nruns = 30
%   log10(I0) transformation applied (same as SM55)
%   Bounds: lb=[0,-12,0,1,1]  ub=[1,-3,1,1000,2]
%
% Easwarakhanthan reference parameters:
%   Iph=0.7608 A, I0=3.223e-7 A, Rs=0.03638 Ohm, Rsh=53.718 Ohm, n=1.4837
%
% Figures saved (for LaTeX manuscript):
%   rtc_iv_fit.eps        — fitted I-V curve, RTH best run
%   rtc_convergence.eps   — mean convergence curves, all algorithms
%   rtc_rmse_boxplot.eps  — RMSE boxplots over 30 runs
% =========================================================================

clear; clc; close all;
rng('default');

%% ===== GLOBAL FIGURE SETTINGS (publication quality) ====================
set(groot, 'DefaultAxesFontSize',    12);
set(groot, 'DefaultAxesFontName',    'Times New Roman');
set(groot, 'DefaultTextFontName',    'Times New Roman');
set(groot, 'DefaultLineLineWidth',    1.8);
set(groot, 'DefaultAxesLineWidth',    1.0);
set(groot, 'DefaultAxesBox',         'on');

%% ===== 1. RTC FRANCE DATASET ============================================
% 26 measured I-V points at T=33 degC
% Source: Easwarakhanthan et al. (1986), Table 1

Ns = 1;
T  = 33 + 273.15;   % [K] — 33 degC as reported in Easwarakhanthan 1986
ref_params_phys = [0.7608, 3.223e-7, 0.03638, 53.718, 1.4837];

% -------------------------------------------------------------------------
% Reference parameters (Easwarakhanthan et al. 1986):
%   Iph=0.7608 A, I0=3.223e-7 A, Rs=0.03638 Ohm, Rsh=53.718 Ohm, n=1.4837
%
% Dataset: 26 measured voltage-current points from -0.2057 V to 0.5900 V.
% Negative-current points near and beyond open circuit are retained because
% they are part of the standard RTC France benchmark table.
% -------------------------------------------------------------------------

% 26 measured voltage-current pairs from the standard RTC France benchmark
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

% Dataset validation
fprintf('=== RTC France 57mm Dataset ===\n');
fprintf('  Points : %d\n', length(V_exp));
fprintf('  V range : [%.4f, %.4f] V\n', V_exp(1), V_exp(end));
fprintf('  Current range : [%.4f, %.4f] A\n', min(I_exp), max(I_exp));
fprintf('  T      = %.2f K (%.1f degC)\n', T, T-273.15);
fprintf('  Reference: Iph=0.7608 A, I0=3.223e-7 A, Rs=0.03638 Ohm, Rsh=53.718 Ohm, n=1.4837\n\n');

%% ===== 2. OPTIMIZATION SETTINGS (same as SM55 benchmark) ================
dim  = 5;
lb   = [0,  -12,   0,     1,   1];   % [Iph, log10(I0), Rs, Rsh, n]
ub   = [1,   -3,   0.5, 500,   2];   % Rsh<500 and Rs<0.5 for single cell
Npop = 50;
Tmax = 1000;
Nruns= 30;

fobj = @(x) PV_objective(x, V_exp, I_exp, Ns, T);

%% ===== 3. ALGORITHM LIST ================================================
algo_names = {'RTH','PSO','GWO','WOA','HHO','AVOA','MGO'};
algo_funcs = {@RTH, @PSO, @GWO, @WOA, @HHO, @AVOA, @MGO};
Nalgo = numel(algo_names);

% Consistent color palette (same as SM55 script)
algo_colors = [
    0.00  0.45  0.74;   % RTH  — blue
    0.85  0.33  0.10;   % PSO  — orange
    0.47  0.67  0.19;   % GWO  — green
    0.49  0.18  0.56;   % WOA  — purple
    0.93  0.69  0.13;   % HHO  — gold
    0.30  0.75  0.93;   % AVOA — cyan
    0.64  0.08  0.18];  % MGO  — dark red

%% ===== 4. MAIN OPTIMIZATION LOOP ========================================
RMSE_all   = zeros(Nalgo, Nruns);
params_all = zeros(Nalgo, Nruns, dim);
conv_all   = zeros(Nalgo, Nruns, Tmax);
cpu_time   = zeros(Nalgo, 1);

fprintf('===========================================================\n');
fprintf(' MULTI-ALGORITHM COMPARISON — RTC France 57mm\n');
fprintf(' %d algorithms x %d runs x Npop=%d x Tmax=%d\n', ...
        Nalgo, Nruns, Npop, Tmax);
fprintf('===========================================================\n\n');

for ia = 1:Nalgo
    fprintf('[%d/%d] %s\n', ia, Nalgo, algo_names{ia});
    t_start = tic;

    for r = 1:Nruns
        rng(ia * 1000 + r);   % same seed scheme as SM55 — reproducible
        [bf, bp, cc] = algo_funcs{ia}(Npop, Tmax, lb, ub, dim, fobj);
        RMSE_all(ia, r)      = bf;
        params_all(ia, r, :) = bp;
        conv_all(ia, r, :)   = cc;
        if mod(r,10)==0
            fprintf('   Run %2d/%d  RMSE = %.4e\n', r, Nruns, bf);
        end
    end

    cpu_time(ia) = toc(t_start);
    fprintf('   CPU: %.1f s  |  RMSE_mean = %.4e\n\n', ...
            cpu_time(ia), mean(RMSE_all(ia,:)));
end

%% ===== 5. PHYSICAL PARAMETERS + METRICS (best run per algorithm) ========
best_params_phys = zeros(Nalgo, dim);
metrics_all      = zeros(Nalgo, 4);     % [RMSE, MAE, MBE, SSE]
best_run_idx     = zeros(Nalgo, 1);

for ia = 1:Nalgo
    [~, br] = min(RMSE_all(ia,:));
    best_run_idx(ia) = br;
    bp = squeeze(params_all(ia, br, :))';
    bp_phys = bp;
    bp_phys(2) = 10^(bp(2));           % log10(I0) -> I0 linear
    best_params_phys(ia,:) = bp_phys;

    I_mod = PV_current_model(bp_phys, V_exp, Ns, T);
    [RMSE, MAE, MBE, SSE] = Compute_metrics(I_exp, I_mod);
    metrics_all(ia,:) = [RMSE, MAE, MBE, SSE];
end

% Statistics over 30 runs
stats = zeros(Nalgo, 4);   % [best, mean, std, worst]
for ia = 1:Nalgo
    rm = RMSE_all(ia,:);
    stats(ia,:) = [min(rm), mean(rm), std(rm), max(rm)];
end

% Average CPU per run
cpu_per_run = cpu_time / Nruns;

%% ===== 6. CONSOLE SUMMARY ===============================================
fprintf('\n============================================================\n');
fprintf(' TABLE — Best Extracted Parameters (Best Run)\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-6s %9s %11s %8s %9s %6s %11s\n', ...
    'Algo','Iph(A)','I0(A)','Rs(Ohm)','Rsh(Ohm)','n','RMSE(A)');
fprintf('------------------------------------------------------------\n');
for ia = 1:Nalgo
    fprintf('%-6s %9.6f %11.4e %8.6f %9.4f %6.4f %11.4e\n', ...
        algo_names{ia}, best_params_phys(ia,1), best_params_phys(ia,2), ...
        best_params_phys(ia,3), best_params_phys(ia,4), ...
        best_params_phys(ia,5), metrics_all(ia,1));
end
fprintf('%-6s %9.6f %11.4e %8.6f %9.4f %6.4f   (Easwarakhanthan ref)\n', ...
    'Ref', 0.7608, 3.223e-7, 0.03638, 53.718, 1.4837);
fprintf('============================================================\n\n');

fprintf('============================================================\n');
fprintf(' TABLE — Statistical Results (30 Independent Runs)\n');
fprintf('------------------------------------------------------------\n');
fprintf('%-6s %11s %11s %11s %11s %8s\n', ...
    'Algo','Best RMSE','Mean RMSE','Std RMSE','Worst RMSE','CPU/run(s)');
fprintf('------------------------------------------------------------\n');
for ia = 1:Nalgo
    fprintf('%-6s %11.4e %11.4e %11.4e %11.4e %8.2f\n', ...
        algo_names{ia}, stats(ia,1), stats(ia,2), ...
        stats(ia,3), stats(ia,4), cpu_per_run(ia));
end
fprintf('============================================================\n\n');

% RTH vs Easwarakhanthan reference
bp_RTH = best_params_phys(1,:);
ref_eas = [0.7608, 3.223e-7, 0.03638, 53.718, 1.4837];
pnames  = {'Iph','I0','Rs','Rsh','n'};
fprintf('--- RTH vs Easwarakhanthan (1986) Reference ---\n');
fprintf('%-6s %12s %12s %10s\n','Param','RTH','Reference','Error%%');
for p=1:5
    err = abs(bp_RTH(p)-ref_eas(p))/ref_eas(p)*100;
    fprintf('%-6s %12.6g %12.6g %10.4f%%\n', pnames{p}, bp_RTH(p), ref_eas(p), err);
end
fprintf('\n');

fprintf('=== COPY THESE VALUES INTO manuscript.tex [TBF] PLACEHOLDERS ===\n');
fprintf('\n--- Tab. rtc_params (Best Extracted Parameters) ---\n');
for ia = 1:Nalgo
    fprintf('  %s: Iph=%.6f  I0=%.4e  Rs=%.6f  Rsh=%.4f  n=%.4f  RMSE=%.4e\n', ...
        algo_names{ia}, best_params_phys(ia,1), best_params_phys(ia,2), ...
        best_params_phys(ia,3), best_params_phys(ia,4), ...
        best_params_phys(ia,5), metrics_all(ia,1));
end
fprintf('\n--- Tab. rtc_stats (30-Run Statistics) ---\n');
for ia = 1:Nalgo
    fprintf('  %s: Best=%.4e  Mean=%.4e  Std=%.4e  Worst=%.4e  CPU=%.2f s\n', ...
        algo_names{ia}, stats(ia,1), stats(ia,2), stats(ia,3), stats(ia,4), ...
        cpu_per_run(ia));
end
fprintf('=================================================================\n\n');

%% ===== 7. FIGURES =======================================================
V_fine = linspace(min(V_exp)-0.05, max(V_exp)+0.02, 300)';

bp_RTH_phys = best_params_phys(1,:);
I_fine_RTH  = PV_current_model(bp_RTH_phys, V_fine, Ns, T);
I_best_RTH  = PV_current_model(bp_RTH_phys, V_exp, Ns, T);

% ----- Figure 1 : I-V fitted curve — RTH (best run) vs Experimental ------
figure(1);
plot(V_exp, I_exp, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', ...
     'DisplayName', 'Benchmark data');
hold on;
plot(V_fine, I_fine_RTH, '-', 'Color', algo_colors(1,:), ...
     'LineWidth', 2.2, 'DisplayName', 'RTH model (best run)');

% Easwarakhanthan reference-parameter curve, shown only for context
I_ref_fine = PV_current_model(ref_params_phys, V_fine, Ns, T);
plot(V_fine, I_ref_fine, 'k--', 'LineWidth', 1.4, ...
     'DisplayName', 'Easwarakhanthan ref. params');

xlabel('Voltage V (V)');
ylabel('Current I (A)');
legend('Location','southwest','FontSize',10,'Box','on');
grid on;
xlim([min(V_exp)-0.05, max(V_exp)+0.01]);

set(gcf,'Position',[50 50 680 440]);
print(gcf,'rtc_iv_fit','-depsc','-r300');
saveas(gcf,'rtc_iv_fit.png');
fprintf('rtc_iv_fit.eps saved.\n');

% ----- Figure 2 : Convergence curves (mean over 30 runs) -----------------
figure(2);
for ia = 1:Nalgo
    conv_mean = mean(squeeze(conv_all(ia,:,:)), 1);
    semilogy(1:Tmax, conv_mean, '-', 'Color', algo_colors(ia,:), ...
             'LineWidth', 1.8, 'DisplayName', algo_names{ia});
    hold on;
end
xlabel('Iteration');
ylabel('Best RMSE (A)');
legend('Location','northeast','FontSize',10,'Box','on');
grid on;
set(gca,'YMinorGrid','on');
title('RTC France 57mm — Mean Convergence over 30 Runs', ...
      'FontSize',11,'FontName','Times New Roman');
set(gcf,'Position',[110 110 700 440]);
print(gcf,'rtc_convergence','-depsc','-r300');
saveas(gcf,'rtc_convergence.png');
fprintf('rtc_convergence.eps saved.\n');

% ----- Figure 3 : RMSE Boxplot (30 runs) ---------------------------------
figure(3);
boxplot(RMSE_all', algo_names, 'Symbol', 'r+', 'Whisker', 1.5);
set(gca,'YScale','log');
yl_min = min(RMSE_all(:)) * 0.5;
yl_max = max(RMSE_all(:)) * 2;
ylim([max(yl_min, 1e-10), yl_max]);
xlabel('Algorithm');
ylabel('RMSE (A)');
grid on;
title('RTC France 57mm — RMSE over 30 Runs', ...
      'FontSize',11,'FontName','Times New Roman');
set(gcf,'Position',[170 170 720 420]);
print(gcf,'rtc_rmse_boxplot','-depsc','-r300');
saveas(gcf,'rtc_rmse_boxplot.png');
fprintf('rtc_rmse_boxplot.eps saved.\n');

fprintf('\n=== Done. Run pdflatex to update the manuscript. ===\n');
fprintf('Move the 3 .eps files into the figures/ subfolder first.\n');

%% ===== EXPORT: LaTeX Table S2 (30-run RMSE for supplementary) =============
fid = fopen('supp_table_S2_rtc.tex', 'w');
fprintf(fid, '%% Table S2 — auto-generated by Main_RTC_France.m\n');
fprintf(fid, '%% Paste the \\midrule ... \\bottomrule block into supplementary.tex\n\n');
fprintf(fid, '\\midrule\n');
for r = 1:Nruns
    fprintf(fid, '%2d', r);
    for ia = 1:Nalgo
        fprintf(fid, ' & %.6e', RMSE_all(ia, r));
    end
    fprintf(fid, ' \\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fclose(fid);
fprintf('supp_table_S2_rtc.tex written — paste into supplementary.tex Table S2.\n');
