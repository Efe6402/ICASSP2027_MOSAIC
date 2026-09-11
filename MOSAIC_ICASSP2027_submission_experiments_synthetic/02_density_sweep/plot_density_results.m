function plot_density_results(T,outdir)
%PLOT_DENSITY_RESULTS Save F1 and AUC component panels versus density.
figdir=fullfile(outdir,'FIGURES'); if ~isfolder(figdir), mkdir(figdir); end
make_one(T,'f1',figdir); make_one(T,'auc',figdir);
end

function make_one(T,metric,figdir)
defs={'supra','Supra';'within','Within-view';'cross','Cross-node/cross-view'; ...
    'copy','Copy';'component','Mean(within,cross,copy)'};
methods=["MOSAIC","PGL2021","ZW2024"];
colors=[0 .447 .741;.850 .325 .098;.929 .694 .125]; styles={'-','--',':'};
f=figure('Color','w','Position',[50 50 1900 1050]);
L=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
title(L,sprintf('p=5000-selected top-3 | fixed validation thresholds | test %s',upper(metric)), ...
    'FontWeight','bold'); handles=gobjects(0); labels=strings(0,1);
for d=1:size(defs,1)
    ax=nexttile(L,d); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    for m=1:numel(methods)
        for rank=1:3
            Q=T(T.method==methods(m)&T.selection_rank==rank,:); if isempty(Q), continue; end
            Q=sortrows(Q,'realized_density');
            mn=sprintf('test_%s_%s_mean',metric,defs{d,1});
            sd=sprintf('test_%s_%s_std',metric,defs{d,1});
            h=errorbar(ax,Q.realized_density,Q.(mn),Q.(sd), ...
                'Color',colors(m,:),'LineStyle',styles{rank},'LineWidth',1.6, ...
                'Marker','o','MarkerFaceColor','w','MarkerSize',5,'CapSize',5);
            if d==1
                handles(end+1)=h; labels(end+1)=sprintf('%s top-%d',methods(m),rank); %#ok<AGROW>
            end
        end
    end
    ylim(ax,[0 1.02]); xlabel(ax,'realized mean modal edge density');
    ylabel(ax,upper(metric)); title(ax,[defs{d,2},' ',upper(metric)]);
end
lgd=legend(handles,cellstr(labels),'FontSize',9); lgd.Layout.Tile=6;
name=['F1_COMPONENT_SELECTED_',upper(metric),'_COMPONENTS_VS_DENSITY'];
exportgraphics(f,fullfile(figdir,[name,'.png']),'Resolution',220);
exportgraphics(f,fullfile(figdir,[name,'.pdf']),'ContentType','vector'); close(f);
end
