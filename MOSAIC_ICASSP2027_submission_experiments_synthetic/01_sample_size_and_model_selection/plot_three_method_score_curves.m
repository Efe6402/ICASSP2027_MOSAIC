function plot_three_method_score_curves(T,outdir)
%PLOT_THREE_METHOD_SCORE_CURVES Test metrics for component-F1 selections.
defs={'test_auc_supra_mean','test_auc_supra_std','Supra AUC','AUC_SUPRA'; ...
    'test_auc_component_mean','test_auc_component_std','Mean component AUC','AUC_COMPONENT'; ...
    'test_f1_supra_mean','test_f1_supra_std','Supra F1','F1_SUPRA'; ...
    'test_f1_component_mean','test_f1_component_std','Mean component F1','F1_COMPONENT'};
methods=["MOSAIC","PGL2021","ZW2024"];
for d=1:size(defs,1)
    f=figure('Visible','off','Color','w','Position',[100 100 1050 650]); ax=axes(f); hold(ax,'on');
    colors=lines(3); styles={'-','--',':'};
    for m=1:3
        for rank=1:3
            Z=T(T.method==methods(m)&T.selection_criterion=="f1_component"&T.selection_rank==rank,:);
            Z=sortrows(Z,'sample_count'); if isempty(Z), continue; end
            errorbar(ax,Z.sample_count,Z.(defs{d,1}),Z.(defs{d,2}), ...
                'Color',colors(m,:),'LineStyle',styles{rank},'Marker','o', ...
                'LineWidth',1.5,'DisplayName',sprintf('%s top-%d',methods(m),rank));
        end
    end
    set(ax,'XScale','log'); grid(ax,'on'); ylim(ax,[0 1.03]);
    xlabel(ax,'number of graph signals p'); ylabel(ax,defs{d,3});
    title(ax,['Component-F1-selected top-3: ',defs{d,3}]); legend(ax,'Location','bestoutside');
    exportgraphics(f,fullfile(outdir,[defs{d,4},'_VS_P_TOP3.png']),'Resolution',200); close(f);
end
end
