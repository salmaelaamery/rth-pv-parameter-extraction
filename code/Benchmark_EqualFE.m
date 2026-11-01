function Benchmark_EqualFE(Nruns, MaxFEs, output_prefix, make_figures)
%BENCHMARK_EQUALFE Equal function-evaluation benchmark for PV extraction.
%
% All algorithms receive exactly 149050 objective-function evaluations:
%   RTH: 50 initial + 1000*(50 + 49 + 50) = 149050
%   Single-evaluation algorithms: 50 initial + 2980*50 = 149050
% HHO has a branch-dependent evaluation rate and therefore uses HHO_FE,
% whose progress variable and stopping condition are both FE-based.

if nargin < 1 || isempty(Nruns), Nruns = 30; end
if nargin < 2 || isempty(MaxFEs), MaxFEs = 149050; end
if nargin < 3 || isempty(output_prefix), output_prefix = ''; end
if nargin < 4 || isempty(make_figures), make_figures = true; end

clc;
if make_figures
    close all;
end

Npop = 50;
dim = 5;

algo_names = {'RTH','PSO','GWO','WOA','HHO','AVOA','MGO'};
algo_funcs = {@RTH,@PSO,@GWO,@WOA,@HHO_FE,@AVOA,@MGO};
rth_iterations = ceil((MaxFEs-Npop)/(3*Npop-1));
standard_iterations = ceil((MaxFEs-Npop)/Npop);
driver_iterations = [rth_iterations, standard_iterations, ...
    standard_iterations, standard_iterations, MaxFEs, ...
    standard_iterations, standard_iterations];

datasets = get_datasets();
Nd = numel(datasets);
Na = numel(algo_names);

candidate_grid = [1, 50, 100:100:1000, 2000:1000:MaxFEs, MaxFEs];
fe_grid = unique(candidate_grid(candidate_grid <= MaxFEs));
Ng = numel(fe_grid);

RMSE_all = nan(Nd, Na, Nruns);
params_all = nan(Nd, Na, Nruns, dim);
cpu_all = nan(Nd, Na, Nruns);
FE_curve_all = nan(Nd, Na, Nruns, Ng);
FE_count_all = zeros(Nd, Na, Nruns);

fprintf('============================================================\n');
fprintf(' EQUAL FUNCTION-EVALUATION BENCHMARK\n');
fprintf(' MaxFEs=%d, Npop=%d, Nruns=%d\n', MaxFEs, Npop, Nruns);
fprintf('============================================================\n\n');

for id = 1:Nd
    data = datasets(id);
    base_fobj = @(x) PV_objective(x, data.V, data.I, data.Ns, data.T);
    fprintf('Dataset %d/%d: %s (%d points)\n', id, Nd, data.name, numel(data.V));

    for ia = 1:Na
        fprintf('  Algorithm %d/%d: %s\n', ia, Na, algo_names{ia});

        for r = 1:Nruns
            % The same run seed is reused across algorithms. After removing
            % the unused initial random matrix in RTH.m, all algorithms start
            % from the same uniformly sampled population for a given dataset.
            rng(50000 + id*1000 + r, 'twister');

            [counted_fobj, get_state] = make_counted_objective(base_fobj, MaxFEs);
            t0 = tic;
            algo_funcs{ia}(Npop, driver_iterations(ia), data.lb, data.ub, dim, counted_fobj);
            elapsed = toc(t0);
            state = get_state();

            if state.fe_count ~= MaxFEs
                error('%s/%s/run %d used %d FEs instead of %d.', ...
                    data.name, algo_names{ia}, r, state.fe_count, MaxFEs);
            end

            RMSE_all(id,ia,r) = state.best_value;
            params_all(id,ia,r,:) = state.best_position;
            cpu_all(id,ia,r) = elapsed;
            FE_count_all(id,ia,r) = state.fe_count;
            FE_curve_all(id,ia,r,:) = sample_best_curve( ...
                state.best_curve, fe_grid);

            if mod(r,10) == 0
                fprintf('    run %02d/%02d: RMSE=%.6e, CPU=%.2f s\n', ...
                    r, Nruns, state.best_value, elapsed);
            end
        end
        fprintf('    mean=%.6e, median=%.6e, CPU/run=%.2f s\n', ...
            mean(squeeze(RMSE_all(id,ia,:))), ...
            median(squeeze(RMSE_all(id,ia,:))), ...
            mean(squeeze(cpu_all(id,ia,:))));
    end
    fprintf('\n');
end

stats = compute_stats(RMSE_all);
paired_stats = compute_paired_stats(RMSE_all, algo_names);
save([output_prefix 'EqualFE_Benchmark_results.mat'], 'datasets', 'algo_names', ...
    'RMSE_all', 'params_all', 'cpu_all', 'FE_curve_all', 'FE_count_all', ...
    'fe_grid', 'stats', 'paired_stats', 'Npop', 'MaxFEs', 'Nruns');

write_run_csv([output_prefix 'EqualFE_run_RMSE.csv'], datasets, algo_names, RMSE_all);
write_summary_csv([output_prefix 'EqualFE_summary.csv'], datasets, algo_names, stats, cpu_all);
write_paired_csv([output_prefix 'EqualFE_paired_tests.csv'], datasets, ...
    algo_names, paired_stats);
write_latex_table([output_prefix 'equal_fe_summary.tex'], datasets, algo_names, stats, cpu_all, ...
    Npop, MaxFEs, Nruns);
if make_figures
    make_convergence_figures(datasets, algo_names, fe_grid, FE_curve_all, ...
        output_prefix, Nruns);
    make_boxplot_figures(datasets, algo_names, RMSE_all, output_prefix, Nruns);
end

fprintf('Saved %sEqualFE_Benchmark_results.mat and equal-FE outputs.\n', ...
    output_prefix);
end

function datasets = get_datasets()
datasets(1).name = 'Shell SM55';
datasets(1).short = 'SM55';
datasets(1).Ns = 36;
datasets(1).T = 25 + 273.15;
datasets(1).V = [0.0000;1.3571;2.7143;4.0714;5.4286;6.7857;8.1429; ...
    9.5000;10.8571;12.2143;13.5714;14.9286;16.2857;17.6429; ...
    19.0000;19.2650;19.5300;19.7950;20.0600;20.3250;20.5900; ...
    20.8550;21.1200;21.3850;21.6500];
datasets(1).I = [3.4500;3.4498;3.4496;3.4493;3.4490;3.4485;3.4478; ...
    3.4463;3.4431;3.4360;3.4198;3.3826;3.2972;3.1029;2.6709; ...
    2.5403;2.3893;2.2152;2.0150;1.7856;1.5239;1.2264;0.8900; ...
    0.5112;0.0871];
datasets(1).lb = [0,-12,0,1,1];
datasets(1).ub = [5,-3,0.5,8000,2];

datasets(2).name = 'RTC France';
datasets(2).short = 'RTC';
datasets(2).Ns = 1;
datasets(2).T = 33 + 273.15;
datasets(2).V = [-0.2057;-0.1291;-0.0588;0.0057;0.0646;0.1185;0.1678; ...
    0.2132;0.2545;0.2924;0.3269;0.3585;0.3873;0.4137;0.4373; ...
    0.4590;0.4784;0.4960;0.5119;0.5265;0.5398;0.5521;0.5633; ...
    0.5736;0.5833;0.5900];
datasets(2).I = [0.7640;0.7620;0.7605;0.7605;0.7600;0.7590;0.7570; ...
    0.7570;0.7555;0.7540;0.7505;0.7465;0.7385;0.7280;0.7065; ...
    0.6755;0.6320;0.5730;0.4990;0.4130;0.3165;0.2120;0.1035; ...
    -0.0100;-0.1230;-0.2100];
datasets(2).lb = [0,-12,0,1,1];
datasets(2).ub = [1,-3,0.5,500,2];

datasets(3).name = 'Photowatt-PWP201';
datasets(3).short = 'PWP201';
datasets(3).Ns = 36;
datasets(3).T = 45 + 273.15;
datasets(3).V = [0.1248;1.8093;3.3511;4.7622;6.0538;7.2364;8.3189; ...
    9.3097;10.2163;11.0449;11.8018;12.4929;13.1231;13.6983; ...
    14.2221;14.6995;15.1346;15.5311;15.8929;16.2229;16.5241; ...
    16.7987;17.0499;17.2793;17.4885];
datasets(3).I = [1.0315;1.0300;1.0260;1.0220;1.0180;1.0155;1.0140; ...
    1.0100;1.0035;0.9880;0.9630;0.9255;0.8725;0.8075;0.7265; ...
    0.6345;0.5345;0.4275;0.3185;0.2085;0.1010;-0.0080;-0.1110; ...
    -0.2090;-0.3030];
datasets(3).lb = [0,-12,0,1,1];
datasets(3).ub = [2,-3,2,5000,2];
end

function sampled = sample_best_curve(curve, fe_grid)
sampled = curve(fe_grid);
end

function stats = compute_stats(RMSE_all)
[Nd,Na,~] = size(RMSE_all);
stats = nan(Nd,Na,6);
for id = 1:Nd
    for ia = 1:Na
        values = squeeze(RMSE_all(id,ia,:));
        stats(id,ia,:) = [min(values), mean(values), median(values), ...
            max(values), std(values), local_iqr(values)];
    end
end
end

function paired_stats = compute_paired_stats(RMSE_all,algo_names)
[Nd,Na,~] = size(RMSE_all);
paired_stats = nan(Nd,Na-1,6);
for id = 1:Nd
    rth = squeeze(RMSE_all(id,1,:));
    raw_p = nan(Na-1,1);
    for ia = 2:Na
        competitor = squeeze(RMSE_all(id,ia,:));
        delta = competitor-rth;
        nonzero = delta(delta ~= 0);
        if isempty(nonzero)
            p_value = 1;
        elseif exist('signrank','file') == 2
            p_value = signrank(competitor,rth,'method','approximate');
        else
            p_value = two_sided_sign_test(nonzero);
        end
        raw_p(ia-1) = p_value;
        paired_stats(id,ia-1,1) = p_value;
        paired_stats(id,ia-1,2) = nan;
        paired_stats(id,ia-1,3) = median(delta);
        paired_stats(id,ia-1,4) = mean(delta > 0);
        paired_stats(id,ia-1,5) = mean(delta < 0);
    end
    adjusted = holm_adjust(raw_p);
    paired_stats(id,:,2) = adjusted';
    paired_stats(id,:,6) = (2:numel(algo_names))';
end
end

function p_value = two_sided_sign_test(nonzero)
n = numel(nonzero);
k = min(sum(nonzero > 0),sum(nonzero < 0));
p_value = min(1,2*sum(arrayfun(@(j) nchoosek(n,j),0:k))/2^n);
end

function adjusted = holm_adjust(raw_p)
[sorted_p,order] = sort(raw_p);
m = numel(raw_p);
sorted_adjusted = zeros(m,1);
running_max = 0;
for i = 1:m
    running_max = max(running_max,(m-i+1)*sorted_p(i));
    sorted_adjusted(i) = min(1,running_max);
end
adjusted = zeros(m,1);
adjusted(order) = sorted_adjusted;
end

function value = local_iqr(values)
values = sort(values(:));
value = local_percentile(values,75) - local_percentile(values,25);
end

function value = local_percentile(values,p)
n = numel(values);
position = 1 + (n-1)*p/100;
lo = floor(position);
hi = ceil(position);
if lo == hi
    value = values(lo);
else
    value = values(lo) + (position-lo)*(values(hi)-values(lo));
end
end

function write_run_csv(filename,datasets,algo_names,RMSE_all)
fid = fopen(filename,'w');
fprintf(fid,'Dataset,Algorithm,Run,RMSE\n');
for id = 1:numel(datasets)
    for ia = 1:numel(algo_names)
        for r = 1:size(RMSE_all,3)
            fprintf(fid,'%s,%s,%d,%.12e\n',datasets(id).short, ...
                algo_names{ia},r,RMSE_all(id,ia,r));
        end
    end
end
fclose(fid);
end

function write_summary_csv(filename,datasets,algo_names,stats,cpu_all)
fid = fopen(filename,'w');
fprintf(fid,'Dataset,Algorithm,Best,Mean,Median,Worst,Std,IQR,MeanCPU_s\n');
for id = 1:numel(datasets)
    for ia = 1:numel(algo_names)
        fprintf(fid,'%s,%s,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.6f\n', ...
            datasets(id).short,algo_names{ia},stats(id,ia,1),stats(id,ia,2), ...
            stats(id,ia,3),stats(id,ia,4),stats(id,ia,5),stats(id,ia,6), ...
            mean(squeeze(cpu_all(id,ia,:))));
    end
end
fclose(fid);
end

function write_paired_csv(filename,datasets,algo_names,paired_stats)
fid = fopen(filename,'w');
fprintf(fid,['Dataset,Comparison,PairedWilcoxonP,HolmAdjustedP,' ...
    'MedianDeltaCompetitorMinusRTH,RTHWinFraction,CompetitorWinFraction\n']);
for id = 1:numel(datasets)
    for ic = 1:(numel(algo_names)-1)
        fprintf(fid,'%s,RTH_vs_%s,%.12e,%.12e,%.12e,%.6f,%.6f\n', ...
            datasets(id).short,algo_names{ic+1},paired_stats(id,ic,1), ...
            paired_stats(id,ic,2),paired_stats(id,ic,3), ...
            paired_stats(id,ic,4),paired_stats(id,ic,5));
    end
end
fclose(fid);
end

function write_latex_table(filename,datasets,algo_names,stats,cpu_all,Npop,MaxFEs,Nruns)
fid = fopen(filename,'w');
fprintf(fid,'%% Auto-generated by Benchmark_EqualFE.m\n');
fprintf(fid,'\\begin{table*}[htbp]\n\\centering\n');
fprintf(fid,'\\caption{Equal function-evaluation comparison over %d runs ($N=%d$, maximum FEs $=%d$).}\n',Nruns,Npop,MaxFEs);
fprintf(fid,'\\label{tab:equal_fe_summary}\n\\small\n');
fprintf(fid,'\\resizebox{\\textwidth}{!}{%%\n');
fprintf(fid,'\\begin{tabular}{@{}llccccccc@{}}\n\\toprule\n');
fprintf(fid,'Benchmark & Algorithm & Best & Mean & Median & Worst & Std. dev. & IQR & CPU/run (s) \\\\\n\\midrule\n');
for id = 1:numel(datasets)
    for ia = 1:numel(algo_names)
        fprintf(fid,'%s & %s & %.4e & %.4e & %.4e & %.4e & %.4e & %.4e & %.2f \\\\\n', ...
            datasets(id).short,algo_names{ia},stats(id,ia,1),stats(id,ia,2), ...
            stats(id,ia,3),stats(id,ia,4),stats(id,ia,5),stats(id,ia,6), ...
            mean(squeeze(cpu_all(id,ia,:))));
    end
    if id < numel(datasets)
        fprintf(fid,'\\midrule\n');
    end
end
fprintf(fid,'\\bottomrule\n\\end{tabular}}\n\\end{table*}\n');
fclose(fid);
end

function make_convergence_figures(datasets,algo_names,fe_grid,FE_curve_all, ...
    output_prefix,Nruns)
colors = lines(numel(algo_names));
for id = 1:numel(datasets)
    fig = figure('Visible','off','Color','w');
    hold on;
    for ia = 1:numel(algo_names)
        curves = squeeze(FE_curve_all(id,ia,:,:));
        median_curve = median(curves,1);
        semilogy(fe_grid,median_curve,'LineWidth',1.7,'Color',colors(ia,:), ...
            'DisplayName',algo_names{ia});
    end
    set(gca,'YScale','log');
    grid on;
    xlabel('Function evaluations');
    ylabel('Median best RMSE (A)');
    title(sprintf('%s: convergence under equal FE budget',datasets(id).name));
    legend('Location','best','NumColumns',2);
    set(gca,'FontName','Times New Roman','FontSize',11);
    set(fig,'Position',[100 100 760 460]);
    subtitle(sprintf('%d independent runs',Nruns));
    print(fig,sprintf('%sEqualFE_Convergence_%s',output_prefix, ...
        datasets(id).short),'-depsc','-r300');
    close(fig);
end
end

function make_boxplot_figures(datasets,algo_names,RMSE_all,output_prefix,Nruns)
for id = 1:numel(datasets)
    fig = figure('Visible','off','Color','w');
    boxplot(squeeze(RMSE_all(id,:,:))',algo_names,'Symbol','r+','Whisker',1.5);
    set(gca,'YScale','log','FontName','Times New Roman','FontSize',11);
    values = squeeze(RMSE_all(id,:,:));
    ylim([0.8*min(values(:)),1.25*max(values(:))]);
    grid on;
    xlabel('Algorithm');
    ylabel('RMSE (A)');
    title(sprintf('%s: RMSE over %d runs at equal FE budget', ...
        datasets(id).name,Nruns));
    set(fig,'Position',[100 100 760 440]);
    print(fig,sprintf('%sEqualFE_Boxplot_%s',output_prefix, ...
        datasets(id).short),'-depsc','-r300');
    close(fig);
end
end
