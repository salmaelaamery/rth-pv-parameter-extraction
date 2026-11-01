function Generate_EqualFE_Supplementary()
%GENERATE_EQUALFE_SUPPLEMENTARY Build the equal-FE supplementary document.

load('EqualFE_Benchmark_results.mat','datasets','algo_names','RMSE_all', ...
    'stats','paired_stats','Npop','MaxFEs','Nruns');

fid = fopen('supplementary.tex','w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid,'\\documentclass[12pt,a4paper]{article}\n');
fprintf(fid,'\\usepackage{booktabs}\n\\usepackage{graphicx}\n');
fprintf(fid,'\\usepackage[margin=2.0cm]{geometry}\n');
fprintf(fid,'\\usepackage[colorlinks=true,linkcolor=blue]{hyperref}\n');
fprintf(fid,'\\begin{document}\n\n');
fprintf(fid,'\\begin{center}\n');
fprintf(fid,'{\\Large\\bfseries Supplementary Material}\\\\[6pt]\n');
fprintf(fid,'{\\large Log-Transformed Red-Tailed Hawk Algorithm for\\\\\n');
fprintf(fid,'Single-Diode Photovoltaic Parameter Extraction}\\\\[4pt]\n');
fprintf(fid,'Salma El-Aamry\\\\\n');
fprintf(fid,['Laboratory of Engineering, Industrial Management and Innovation,\\\\\n' ...
    'Faculty of Sciences and Technology, Sidi Mohamed Ben Abdellah University,\\\\\n' ...
    'Fez, Morocco\n']);
fprintf(fid,'\\end{center}\n\n');

fprintf(fid,['This document reports the complete equal function-evaluation (FE) results ' ...
    'used in the revised manuscript.  Each algorithm receives exactly ' ...
    '$%d$ objective evaluations per run, with $N=%d$ and %d matched runs. ' ...
    'For dataset $d$ and run $r$, all algorithms use ' ...
    '\\texttt{rng(50000 + 1000$d$ + $r$)}, so the initial population is ' ...
    'identical within a matched run.  RTH uses 1000 iterations; the ' ...
    'single-evaluation algorithms use 2980 iterations; HHO uses an FE-based ' ...
    'progress and stopping rule because its evaluation count per iteration ' ...
    'is branch-dependent.\n\n'],MaxFEs,Npop,Nruns);
fprintf(fid,['The near-zero RTH dispersion is reported as numerical convergence ' ...
    'under this protocol, not as proof of a certified global optimum.  The ' ...
    'MATLAB files \\texttt{Benchmark\\_EqualFE.m}, ' ...
    '\\texttt{make\\_counted\\_objective.m}, \\texttt{HHO\\_FE.m}, and the ' ...
    'saved MAT/CSV outputs accompany this document.\n\n']);

write_summary_table(fid,datasets,algo_names,stats);
for id = 1:numel(datasets)
    write_runs_table(fid,datasets(id),algo_names,squeeze(RMSE_all(id,:,:)),id);
end
write_paired_table(fid,datasets,algo_names,paired_stats);
write_checklist(fid);

fprintf(fid,'\\end{document}\n');
clear cleanup;
end

function write_summary_table(fid,datasets,algo_names,stats)
fprintf(fid,'\\section*{Table S1 -- Equal-FE statistical summary}\n');
fprintf(fid,'\\begin{table}[htbp]\\centering\\small\n');
fprintf(fid,['\\caption{RMSE summary over 30 matched runs.  All RMSE values are ' ...
    'in amperes.}\\label{tab:supp_summary}\n']);
fprintf(fid,'\\resizebox{\\textwidth}{!}{%%\n');
fprintf(fid,'\\begin{tabular}{llcccccc}\\toprule\n');
fprintf(fid,'Dataset & Algorithm & Best & Mean & Median & Worst & Std. dev. & IQR \\\\\n');
fprintf(fid,'\\midrule\n');
for id = 1:numel(datasets)
    for ia = 1:numel(algo_names)
        fprintf(fid,'%s & %s & %.4e & %.4e & %.4e & %.4e & %.4e & %.4e \\\\\n', ...
            datasets(id).short,algo_names{ia},stats(id,ia,1),stats(id,ia,2), ...
            stats(id,ia,3),stats(id,ia,4),stats(id,ia,5),stats(id,ia,6));
    end
    if id < numel(datasets), fprintf(fid,'\\midrule\n'); end
end
fprintf(fid,'\\bottomrule\\end{tabular}}\n\\end{table}\n\\clearpage\n\n');
end

function write_runs_table(fid,dataset,algo_names,values,table_index)
fprintf(fid,'\\section*{Table S%d -- Complete RMSE values: %s}\n', ...
    table_index+1,dataset.name);
fprintf(fid,'\\begin{table}[htbp]\\centering\\small\n');
fprintf(fid,['\\caption{Complete equal-FE RMSE values (A) for %s over 30 matched ' ...
    'runs.}\\label{tab:supp_%s_runs}\n'],dataset.name,lower(dataset.short));
fprintf(fid,'\\resizebox{\\textwidth}{!}{%%\n');
fprintf(fid,'\\begin{tabular}{cccccccc}\\toprule\n');
fprintf(fid,'Run');
for ia = 1:numel(algo_names), fprintf(fid,' & %s',algo_names{ia}); end
fprintf(fid,' \\\\\n\\midrule\n');
for r = 1:size(values,2)
    fprintf(fid,'%d',r);
    for ia = 1:size(values,1)
        fprintf(fid,' & %.8e',values(ia,r));
    end
    fprintf(fid,' \\\\\n');
end
fprintf(fid,'\\bottomrule\\end{tabular}}\n\\end{table}\n\\clearpage\n\n');
end

function write_paired_table(fid,datasets,algo_names,paired_stats)
fprintf(fid,'\\section*{Table S5 -- Paired statistical comparisons}\n');
fprintf(fid,'\\begin{table}[htbp]\\centering\\small\n');
fprintf(fid,['\\caption{Paired Wilcoxon signed-rank comparisons of RTH with each ' ...
    'competitor under matched initial populations.}\\label{tab:supp_paired}\n']);
fprintf(fid,'\\resizebox{\\textwidth}{!}{%%\n');
fprintf(fid,'\\begin{tabular}{llccccc}\\toprule\n');
fprintf(fid,['Dataset & Comparison & $p$ & Holm $p$ & Median $\\Delta$ & ' ...
    'RTH wins & Competitor wins \\\\\n\\midrule\n']);
for id = 1:numel(datasets)
    for ic = 1:(numel(algo_names)-1)
        fprintf(fid,'%s & RTH vs. %s & %.3e & %.3e & %.4e & %.2f & %.2f \\\\\n', ...
            datasets(id).short,algo_names{ic+1},paired_stats(id,ic,1), ...
            paired_stats(id,ic,2),paired_stats(id,ic,3), ...
            paired_stats(id,ic,4),paired_stats(id,ic,5));
    end
    if id < numel(datasets), fprintf(fid,'\\midrule\n'); end
end
fprintf(fid,'\\bottomrule\\end{tabular}}\n');
fprintf(fid,['\\par\\medskip\\footnotesize $\\Delta$ is competitor RMSE minus RTH ' ...
    'RMSE.  The signed-rank values are supporting evidence because the RTH ' ...
    'samples have negligible dispersion at the displayed precision.\n']);
fprintf(fid,'\\end{table}\n\\clearpage\n\n');
end

function write_checklist(fid)
fprintf(fid,'\\section*{Supplementary reproducibility checklist}\n');
fprintf(fid,'\\begin{itemize}\n');
fprintf(fid,['\\item Complete run-level values are provided above and in ' ...
    '\\texttt{EqualFE\\_run\\_RMSE.csv}.\n']);
fprintf(fid,['\\item Summary statistics and paired tests are provided in ' ...
    '\\texttt{EqualFE\\_summary.csv} and ' ...
    '\\texttt{EqualFE\\_paired\\_tests.csv}.\n']);
fprintf(fid,['\\item Best physical parameters and error metrics are provided in ' ...
    '\\texttt{EqualFE\\_best\\_parameters.csv}.\n']);
fprintf(fid,['\\item The complete best-so-far traces are stored in ' ...
    '\\texttt{EqualFE\\_Benchmark\\_results.mat}.\n']);
fprintf(fid,['\\item The three benchmark datasets, parameter bounds, FE counter, and ' ...
    'seed protocol are contained in \\texttt{Benchmark\\_EqualFE.m}.\n']);
fprintf(fid,['\\item The log-versus-linear $I_0$ ablation code and outputs are ' ...
    'included separately.\n']);
fprintf(fid,['\\item A permanent public repository URL and software license should be ' ...
    'added before final publication.\n']);
fprintf(fid,'\\end{itemize}\n\n');
end
