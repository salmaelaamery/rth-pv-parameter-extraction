% =========================================================================
% Ablation_LogI0_RTH.m
% =========================================================================
% RTH ablation for the saturation-current search variable.
%
% Compared variants:
%   1) log10(I0): x(2) = log10(I0), evaluated with PV_objective.m
%   2) linear I0: x(2) = I0, evaluated with the local objective below
%
% Both variants use the same RTH implementation, Newton-Raphson current solver,
% population size, iteration count, and run seeds.
% =========================================================================

clear; clc; close all;

Npop  = 50;
Tmax  = 1000;
Nruns = 30;
dim   = 5;

datasets = {'SM55','RTC France'};
Nd = numel(datasets);
Nm = 2;
method_names = {'log10(I0)','linear I0'};

RMSE_all = zeros(Nd, Nm, Nruns);
BestPos_all = zeros(Nd, Nm, Nruns, dim);
Conv_all = zeros(Nd, Nm, Nruns, Tmax);
CPU_all = zeros(Nd, Nm);

fprintf('============================================================\n');
fprintf(' RTH ABLATION: log10(I0) vs linear I0\n');
fprintf(' Npop=%d, Tmax=%d, Nruns=%d\n', Npop, Tmax, Nruns);
fprintf('============================================================\n\n');

for id = 1:Nd
    [V_exp, I_exp, Ns, T, lb_log, ub_log, lb_lin, ub_lin] = get_dataset(datasets{id});

    fprintf('Dataset: %s  (%d points, Ns=%d, T=%.2f K)\n', ...
        datasets{id}, numel(V_exp), Ns, T);

    fobj_log = @(x) PV_objective(x, V_exp, I_exp, Ns, T);
    fobj_lin = @(x) PV_objective_linearI0(x, V_exp, I_exp, Ns, T);

    for im = 1:Nm
        if im == 1
            lb = lb_log; ub = ub_log; fobj = fobj_log;
        else
            lb = lb_lin; ub = ub_lin; fobj = fobj_lin;
        end

        fprintf('  Method: %s\n', method_names{im});
        t0 = tic;
        for r = 1:Nruns
            rng(9000 + id*1000 + r);  % same seed sequence for log and linear variants
            [bf, bp, cc] = RTH(Npop, Tmax, lb, ub, dim, fobj);
            RMSE_all(id, im, r) = bf;
            BestPos_all(id, im, r, :) = bp;
            Conv_all(id, im, r, :) = cc;
            if mod(r, 10) == 0
                fprintf('    run %02d/%02d: RMSE = %.6e\n', r, Nruns, bf);
            end
        end
        CPU_all(id, im) = toc(t0);
        fprintf('    mean RMSE = %.6e, CPU/run = %.2f s\n\n', ...
            mean(RMSE_all(id, im, :)), CPU_all(id, im)/Nruns);
    end
end

%% Summary statistics
stats = zeros(Nd, Nm, 6);  % best, mean, median, worst, std, IQR
for id = 1:Nd
    for im = 1:Nm
        rm = squeeze(RMSE_all(id, im, :));
        stats(id, im, :) = [min(rm), mean(rm), median(rm), max(rm), std(rm), iqr_local(rm)];
    end
end

fprintf('\n============================================================\n');
fprintf(' ABLATION SUMMARY\n');
fprintf('============================================================\n');
for id = 1:Nd
    fprintf('\n%s\n', datasets{id});
    fprintf('%-12s %12s %12s %12s %12s %12s %12s %10s\n', ...
        'Method','Best','Mean','Median','Worst','Std','IQR','CPU/run');
    for im = 1:Nm
        fprintf('%-12s %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %10.2f\n', ...
            method_names{im}, stats(id,im,1), stats(id,im,2), stats(id,im,3), ...
            stats(id,im,4), stats(id,im,5), stats(id,im,6), CPU_all(id,im)/Nruns);
    end
end

save('Ablation_LogI0_RTH_results.mat', 'datasets', 'method_names', 'RMSE_all', ...
     'BestPos_all', 'Conv_all', 'CPU_all', 'stats', 'Npop', 'Tmax', 'Nruns');

write_latex_table('ablation_log_i0_results.tex', datasets, method_names, stats, CPU_all, Nruns);
fprintf('\nSaved: Ablation_LogI0_RTH_results.mat\n');
fprintf('Saved: ablation_log_i0_results.tex\n');

%% Local functions
function [V_exp, I_exp, Ns, T, lb_log, ub_log, lb_lin, ub_lin] = get_dataset(name)
    switch name
        case 'SM55'
            Ns = 36;
            T  = 25 + 273.15;
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
            lb_log = [0, -12, 0,   1,    1];
            ub_log = [5,  -3, 0.5, 8000, 2];
            lb_lin = [0, 1e-12, 0,   1,    1];
            ub_lin = [5, 1e-3,  0.5, 8000, 2];

        case 'RTC France'
            Ns = 1;
            T  = 33 + 273.15;
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
            lb_log = [0, -12, 0,   1,   1];
            ub_log = [1,  -3, 0.5, 500, 2];
            lb_lin = [0, 1e-12, 0,   1,   1];
            ub_lin = [1, 1e-3,  0.5, 500, 2];
    end
end

function fitness = PV_objective_linearI0(x, V_exp, I_exp, Ns, T)
    x_phys = x;
    if x_phys(1) < 0 || x_phys(2) <= 0 || x_phys(3) < 0 || x_phys(4) <= 0 || x_phys(5) <= 0
        fitness = 1e10;
        return;
    end

    try
        I_model = PV_current_model(x_phys, V_exp, Ns, T);
    catch
        fitness = 1e10;
        return;
    end

    if any(~isfinite(I_model))
        fitness = 1e10;
        return;
    end

    residuals = I_exp(:) - I_model(:);
    fitness = sqrt(sum(residuals.^2) / numel(I_exp));
    if ~isfinite(fitness)
        fitness = 1e10;
    end
end

function y = iqr_local(x)
    x = sort(x(:));
    y = prctile_local(x, 75) - prctile_local(x, 25);
end

function q = prctile_local(x, p)
    n = numel(x);
    if n == 1
        q = x;
        return;
    end
    pos = 1 + (n - 1) * p / 100;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q = x(lo);
    else
        q = x(lo) + (pos - lo) * (x(hi) - x(lo));
    end
end

function write_latex_table(fname, datasets, method_names, stats, CPU_all, Nruns)
    latex_method_names = {'$\log_{10}(I_0)$', 'linear $I_0$'};
    fid = fopen(fname, 'w');
    fprintf(fid, '%% Auto-generated by Ablation_LogI0_RTH.m\n');
    fprintf(fid, '\\begin{table}[H]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{Ablation of the saturation-current search variable for RTH (%d runs, $N=50$, $T_{\\max}=1000$).}\n', Nruns);
    fprintf(fid, '\\label{tab:ablation_log_i0}\n');
    fprintf(fid, '\\small\n');
    fprintf(fid, '\\resizebox{\\linewidth}{!}{%%\n');
    fprintf(fid, '\\begin{tabular}{@{}llccccccc@{}}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Benchmark & Search variable & Best & Mean & Median & Worst & Std. dev. & IQR & CPU/run (s) \\\\\n');
    fprintf(fid, '\\midrule\n');
    for id = 1:numel(datasets)
        for im = 1:numel(method_names)
            fprintf(fid, '%s & %s & %.4e & %.4e & %.4e & %.4e & %.4e & %.4e & %.2f \\\\\n', ...
                datasets{id}, latex_method_names{im}, stats(id,im,1), stats(id,im,2), ...
                stats(id,im,3), stats(id,im,4), stats(id,im,5), stats(id,im,6), ...
                CPU_all(id,im)/Nruns);
        end
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);
end
