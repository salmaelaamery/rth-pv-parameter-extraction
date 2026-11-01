function Postprocess_EqualFE()
%POSTPROCESS_EQUALFE Regenerate equal-FE figures and export physical results.

load('EqualFE_Benchmark_results.mat','datasets','algo_names','RMSE_all', ...
    'params_all','FE_curve_all','fe_grid','Nruns','MaxFEs');

colors = lines(numel(algo_names));
for id = 1:numel(datasets)
    fig = figure('Visible','off','Color','w');
    hold on;
    for ia = 1:numel(algo_names)
        curves = squeeze(FE_curve_all(id,ia,:,:));
        median_curve = median(curves,1);
        plot(fe_grid,median_curve,'LineWidth',1.7,'Color',colors(ia,:), ...
            'DisplayName',algo_names{ia});
    end
    set(gca,'YScale','log','FontName','Times New Roman','FontSize',11);
    grid on;
    xlabel('Function evaluations');
    ylabel('Median best RMSE (A)');
    title(sprintf('%s: convergence under equal FE budget',datasets(id).name));
    subtitle(sprintf('%d independent runs; maximum FEs = %d',Nruns,MaxFEs));
    legend('Location','best','NumColumns',2);
    set(fig,'Position',[100 100 760 460]);
    print(fig,sprintf('EqualFE_Convergence_%s',datasets(id).short), ...
        '-depsc','-r300');
    close(fig);

    fig = figure('Visible','off','Color','w');
    values = squeeze(RMSE_all(id,:,:));
    boxplot(values',algo_names,'Symbol','r+','Whisker',1.5);
    set(gca,'YScale','log','FontName','Times New Roman','FontSize',11);
    ylim([0.8*min(values(:)),1.25*max(values(:))]);
    grid on;
    xlabel('Algorithm');
    ylabel('RMSE (A)');
    title(sprintf('%s: RMSE over %d runs at equal FE budget', ...
        datasets(id).name,Nruns));
    set(fig,'Position',[100 100 760 440]);
    print(fig,sprintf('EqualFE_Boxplot_%s',datasets(id).short), ...
        '-depsc','-r300');
    close(fig);
end

export_parameters(datasets,algo_names,RMSE_all,params_all);
make_pwp201_fit(datasets,RMSE_all,params_all);
end

function export_parameters(datasets,algo_names,RMSE_all,params_all)
fid = fopen('EqualFE_best_parameters.csv','w');
fprintf(fid,['Dataset,Algorithm,Run,Iph_A,I0_A,Rs_Ohm,Rsh_Ohm,n,' ...
    'RMSE_A,MAE_A,MBE_A,SSE_A2\n']);
for id = 1:numel(datasets)
    for ia = 1:numel(algo_names)
        values = squeeze(RMSE_all(id,ia,:));
        [~,run_index] = min(values);
        x = squeeze(params_all(id,ia,run_index,:))';
        x(2) = 10^x(2);
        modeled = PV_current_model(x,datasets(id).V,datasets(id).Ns, ...
            datasets(id).T);
        [rmse,mae,mbe,sse] = Compute_metrics(datasets(id).I,modeled);
        fprintf(fid,['%s,%s,%d,%.12e,%.12e,%.12e,%.12e,%.12e,' ...
            '%.12e,%.12e,%.12e,%.12e\n'],datasets(id).short, ...
            algo_names{ia},run_index,x(1),x(2),x(3),x(4),x(5), ...
            rmse,mae,mbe,sse);
    end
end
fclose(fid);
end

function make_pwp201_fit(datasets,RMSE_all,params_all)
id = find(strcmp({datasets.short},'PWP201'),1);
[~,run_index] = min(squeeze(RMSE_all(id,1,:)));
x = squeeze(params_all(id,1,run_index,:))';
x(2) = 10^x(2);

voltage_dense = linspace(min(datasets(id).V),max(datasets(id).V),400)';
current_dense = PV_current_model(x,voltage_dense,datasets(id).Ns, ...
    datasets(id).T);
current_points = PV_current_model(x,datasets(id).V,datasets(id).Ns, ...
    datasets(id).T);

fig = figure('Visible','off','Color','w');
subplot(1,2,1);
plot(voltage_dense,current_dense,'b-','LineWidth',1.8);
hold on;
plot(datasets(id).V,datasets(id).I,'ro','MarkerSize',5, ...
    'LineWidth',1.1);
grid on;
xlabel('Voltage (V)');
ylabel('Current (A)');
legend('RTH-SDM','Measured','Location','southwest');
title('(a) Current-voltage fit');

subplot(1,2,2);
plot(voltage_dense,voltage_dense.*current_dense,'b-','LineWidth',1.8);
hold on;
plot(datasets(id).V,datasets(id).V.*datasets(id).I,'ro', ...
    'MarkerSize',5,'LineWidth',1.1);
grid on;
xlabel('Voltage (V)');
ylabel('Power (W)');
legend('RTH-SDM','Measured','Location','southwest');
title('(b) Power-voltage fit');

set(findall(fig,'-property','FontName'),'FontName','Times New Roman');
set(findall(fig,'-property','FontSize'),'FontSize',11);
set(fig,'Position',[100 100 920 400]);
print(fig,'EqualFE_PWP201_Fit','-depsc','-r300');
close(fig);
end
