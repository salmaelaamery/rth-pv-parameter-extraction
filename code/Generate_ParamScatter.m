function Generate_ParamScatter()
%GENERATE_PARAMSCATTER Scatter plots of extracted parameters over 30 runs.

load('EqualFE_Benchmark_results.mat','datasets','algo_names','params_all');

param_labels = {'I_{ph} (A)','log_{10}(I_0)','R_s (\Omega)', ...
    'R_{sh} (\Omega)','n'};
colors = lines(numel(algo_names));

for id = 1:numel(datasets)
    fig = figure('Visible','off','Color','w');
    for ip = 1:5
        subplot(2,3,ip);
        hold on;
        for ia = 1:numel(algo_names)
            values = squeeze(params_all(id,ia,:,ip));
            x = ia + 0.18*(rand(size(values))-0.5);
            scatter(x,values,18,'MarkerFaceColor',colors(ia,:), ...
                'MarkerEdgeColor','k','MarkerFaceAlpha',0.65, ...
                'MarkerEdgeAlpha',0.35);
        end
        grid on;
        xlim([0.4,numel(algo_names)+0.6]);
        set(gca,'XTick',1:numel(algo_names),'XTickLabel',algo_names, ...
            'FontName','Times New Roman','FontSize',9);
        xtickangle(35);
        ylabel(param_labels{ip});
        if ip == 4
            set(gca,'YScale','log');
        end
    end
    sgtitle(sprintf('%s: extracted parameter distributions over 30 runs', ...
        datasets(id).name),'FontName','Times New Roman','FontWeight','bold');
    set(fig,'Position',[100 100 980 560]);
    print(fig,sprintf('figures/EqualFE_ParamScatter_%s',datasets(id).short), ...
        '-depsc','-r300');
    close(fig);
end
end
