function plot_noise_results(S,cfg,out)
%PLOT_NOISE_RESULTS Visualize the two-dimensional kappa x SNR experiment.
if isempty(S), return; end
figdir=fullfile(out,'FIGURES'); if ~isfolder(figdir), mkdir(figdir); end
for family={ 'auc','f1' }
    fam=family{1};
    for s=cfg.plot_snr_db
        if any(abs(cfg.snr_db(isfinite(cfg.snr_db))-s)<1e-12)
            make_kappa_panel(S,fam,figdir,cfg,s);
        end
    end
    make_clean_check(S,fam,figdir,cfg);
    for method=["MOSAIC","PGL2021","ZW2024"]
        make_top1_heatmap(S,fam,figdir,cfg,method);
    end
end
end

function make_kappa_panel(S,family,out,cfg,snr_value)
suffix={'supra','within','cross','copy','component'};
titles={'Supra','Within-view','Cross-node/cross-view','Copy','Mean(within,cross,copy)'};
f=figure('Visible','off','Color','w','Position',[20 20 1850 1050]);
tl=tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
title(tl,sprintf(['p=%d | fixed p=5000 top-%d configurations | SNR=%g dB | ', ...
    'view-wise heterogeneity sweep | %s'],cfg.anchor_p,cfg.top_k,snr_value,upper(family)), ...
    'FontWeight','bold');
methods=["MOSAIC","PGL2021","ZW2024"];
colors=[0 .447 .741;.85 .325 .098;.929 .694 .125];
styles={'-','--',':'};
for k=1:5
    ax=nexttile; hold(ax,'on'); grid(ax,'on');
    for mi=1:3
        for rank=1:cfg.top_k
            Q=S(S.method==methods(mi)&S.selection_rank==rank& ...
                abs(S.target_snr_db-snr_value)<1e-12,:);
            if isempty(Q), continue; end
            Q=sortrows(Q,'kappa');
            meanvar=sprintf('%s_%s_mean',family,suffix{k});
            stdvar=sprintf('%s_%s_std',family,suffix{k});
            errorbar(ax,Q.kappa,Q.(meanvar),Q.(stdvar), ...
                'Color',colors(mi,:),'LineStyle',styles{rank}, ...
                'Marker','o','LineWidth',1.4,'MarkerSize',5);
        end
    end
    title(ax,[titles{k},' ',upper(family)]);
    xlabel(ax,'view heterogeneity \kappa (0 = IID, 1 = endpoint)');
    ylabel(ax,upper(family));
    xlim(ax,[0 1]); xticks(ax,cfg.kappa_grid); ylim(ax,[0 1.02]);
end
ax=nexttile; hold(ax,'on'); handles=gobjects(0); labels=strings(0);
for mi=1:3
    for rank=1:cfg.top_k
        handles(end+1)=plot(ax,nan,nan,'Color',colors(mi,:), ...
            'LineStyle',styles{rank},'Marker','o','LineWidth',1.4); %#ok<AGROW>
        labels(end+1)=methods(mi)+" top-"+rank; %#ok<AGROW>
    end
end
axis(ax,'off'); legend(ax,handles,labels,'Location','northwest');
name=sprintf('FIXED_P5000_TOP%d_%s_COMPONENTS_VS_KAPPA_AT_SNR_%02dDB', ...
    cfg.top_k,upper(family),round(snr_value));
exportgraphics(f,fullfile(out,[name,'.png']),'Resolution',220);
exportgraphics(f,fullfile(out,[name,'.pdf']),'ContentType','vector'); close(f);
end

function make_top1_heatmap(S,family,out,cfg,method)
Q=S(S.method==method&S.selection_rank==1,:);
M=nan(numel(cfg.kappa_grid),numel(cfg.snr_db));
metric=[family,'_supra_mean'];
for ki=1:numel(cfg.kappa_grid)
    for si=1:numel(cfg.snr_db)
        if isinf(cfg.snr_db(si))
            Z=Q(abs(Q.kappa-cfg.kappa_grid(ki))<1e-12&isinf(Q.target_snr_db),:);
        else
            Z=Q(abs(Q.kappa-cfg.kappa_grid(ki))<1e-12& ...
                abs(Q.target_snr_db-cfg.snr_db(si))<1e-12,:);
        end
        if ~isempty(Z)
            vals=Z.(metric); M(ki,si)=vals(1);
        end
    end
end
f=figure('Visible','off','Color','w','Position',[50 50 1050 620]);
ax=axes(f); imagesc(ax,M,[0 1]); colorbar(ax);
xticks(ax,1:numel(cfg.snr_db)); labels=string(cfg.snr_db); labels(1)="clean"; xticklabels(ax,labels);
yticks(ax,1:numel(cfg.kappa_grid)); yticklabels(ax,string(cfg.kappa_grid));
xlabel(ax,'SNR (dB)'); ylabel(ax,'view heterogeneity \kappa');
title(ax,sprintf('%s top-1 | supra %s | view-wise kappa x SNR',method,upper(family)));
name=sprintf('%s_TOP1_SUPRA_%s_KAPPA_X_SNR_HEATMAP',method,upper(family));
exportgraphics(f,fullfile(out,[name,'.png']),'Resolution',220);
exportgraphics(f,fullfile(out,[name,'.pdf']),'ContentType','vector'); close(f);
end

function make_clean_check(S,family,out,cfg)
% Clean endpoint must be invariant to kappa because alpha=0.
Q=S(S.selection_rank==1&isinf(S.target_snr_db),:);
methods=["MOSAIC","PGL2021","ZW2024"];
colors=[0 .447 .741;.85 .325 .098;.929 .694 .125];
f=figure('Visible','off','Color','w','Position',[50 50 900 620]);
ax=axes(f); hold(ax,'on'); grid(ax,'on'); metric=[family,'_supra_mean'];
for mi=1:3
    Z=sortrows(Q(Q.method==methods(mi),:),'kappa');
    plot(ax,Z.kappa,Z.(metric),'Color',colors(mi,:),'Marker','o','LineWidth',1.6);
end
xlabel(ax,'\kappa'); ylabel(ax,['Clean supra ',upper(family)]); ylim(ax,[0 1.02]);
xlim(ax,[0 1]); xticks(ax,cfg.kappa_grid);
title(ax,'Clean endpoint sanity check: performance should not depend on \kappa');
legend(ax,methods,'Location','best');
name=sprintf('CLEAN_ENDPOINT_TOP1_SUPRA_%s_VS_KAPPA_SANITY_CHECK',upper(family));
exportgraphics(f,fullfile(out,[name,'.png']),'Resolution',220);
exportgraphics(f,fullfile(out,[name,'.pdf']),'ContentType','vector'); close(f);
end
