function plot_rho_results(S,cfg,out)
%PLOT_RHO_RESULTS Component AUC and F1 panels for each selection criterion.
if isempty(S), return; end
figdir=fullfile(out,'FIGURES'); if ~isfolder(figdir), mkdir(figdir); end
for criterion=string(cfg.selection_criteria)
    Q=S(S.selection_criterion==criterion,:);
    if isempty(Q), continue; end
    make_panel(Q,'auc',criterion,figdir,cfg.anchor_p,cfg.top_k);
    make_panel(Q,'f1',criterion,figdir,cfg.anchor_p,cfg.top_k);
end
end

function make_panel(S,family,criterion,out,p,top_k)
suffix={'supra','within','cross','copy','component'};
titles={'Supra','Within-view','Cross-node/cross-view','Copy', ...
    'Mean(within,cross,copy)'};
f=figure('Visible','off','Color','w','Position',[20 20 1850 1050]);
tl=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('p=%d | %s-selected top-3 | test %s vs \\rho', ...
    p,strrep(char(criterion),'_',' '),upper(family)),'FontWeight','bold');
methods=["MOSAIC","PGL2021","ZW2024"];
colors=[0 .447 .741;.85 .325 .098;.929 .694 .125];
styles={'-','--',':'};
for k=1:5
    ax=nexttile; hold(ax,'on'); grid(ax,'on');
    for mi=1:3
        for rank=1:top_k
            Q=S(S.method==methods(mi)&S.selection_rank==rank,:);
            if isempty(Q), continue; end
            Q=sortrows(Q,'rho');
            meanvar=sprintf('%s_%s_mean',family,suffix{k});
            stdvar=sprintf('%s_%s_std',family,suffix{k});
            errorbar(ax,Q.rho,Q.(meanvar),Q.(stdvar), ...
                'Color',colors(mi,:),'LineStyle',styles{rank}, ...
                'Marker','o','LineWidth',1.4,'MarkerSize',5);
        end
    end
    title(ax,[titles{k},' ',upper(family)]); xlabel(ax,'profile mixing \rho');
    ylabel(ax,upper(family)); ylim(ax,[0 1.02]); xlim(ax,[min(S.rho) max(S.rho)]);
end
ax=nexttile; hold(ax,'on'); handles=gobjects(0); labels=strings(0);
for mi=1:3
    for rank=1:top_k
        handles(end+1)=plot(ax,nan,nan,'Color',colors(mi,:), ...
            'LineStyle',styles{rank},'Marker','o','LineWidth',1.4); %#ok<AGROW>
        labels(end+1)=methods(mi)+" top-"+rank; %#ok<AGROW>
    end
end
axis(ax,'off'); legend(ax,handles,labels,'Location','northwest');
name=sprintf('%s_SELECTED_%s_COMPONENTS_VS_RHO',upper(char(criterion)),upper(family));
exportgraphics(f,fullfile(out,[name,'.png']),'Resolution',220);
exportgraphics(f,fullfile(out,[name,'.pdf']),'ContentType','vector'); close(f);
end
