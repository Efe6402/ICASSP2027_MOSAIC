function save_selected_recovery_figures(method,S,bank,truth,cfg,p,pdir)
%SAVE_SELECTED_RECOVERY_FIGURES Save all selected fits; plot a clear subset.
bundle=fullfile(pdir,'selected_bundles',char(method));
if ~exist(bundle,'dir'), mkdir(bundle); end
seen=[];
for q=1:height(S)
    cid=double(S.candidate_id(q));
    cache=fullfile(pdir,'fit_cache',char(method),'test', ...
        sprintf('candidate_%04d_trial_%02d.mat',cid,double(bank.trial(1))));
    if ~isfile(cache), continue; end
    e=load(cache,'est'); est=e.est;
    tau=struct('supra',S.tau_supra(q),'within',S.tau_within(q), ...
        'cross',S.tau_cross(q),'copy',S.tau_copy(q));
    out=fullfile(bundle,sprintf('%s_rank%d_candidate%04d', ...
        char(S.selection_criterion(q)),double(S.selection_rank(q)),cid));
    if cfg.save_all_selected_fits
        save([out,'.mat'],'est','tau','-v7.3');
    end
    if S.selection_rank(q)~=1 || ...
            ~ismember(string(S.selection_criterion(q)),["auc_supra","f1_component"])
        continue
    end
    reps=unique([cfg.sample_counts(1),1000,cfg.sample_counts(end)]);
    if ~ismember(p,reps), continue; end
    key=cid+1e6*double(S.selection_criterion(q)=="f1_component");
    if ismember(key,seen), continue; end; seen(end+1)=key; %#ok<AGROW>
    plot_components(est,truth,tau,method,S.selection_criterion(q),out);
    if method=="MOSAIC", plot_mosaic_modes(est,truth,out); end
    if method=="ZW2024", plot_zw_primitive(est,truth,out); end
end
end

function plot_components(est,t,tau,method,criterion,out)
defs={'Full supra','all','all','supra';'Within-view','within','within','within'; ...
    'Cross-node/cross-view','crossnode_crossview','crossnode_crossview','cross'; ...
    'Copy','copy','copy','copy'};
A=est.A_primary; f=figure('Visible','off','Color','w','Position',[50 50 1500 720]);
tl=tiledlayout(f,2,4,'TileSpacing','compact','Padding','compact');
for j=1:4
    mask=t.truth_masks.(defs{j,3}); T=t.A_supra_true.*mask; E=A.*mask;
    lim=max([T(:);E(:)]); if lim<=0, lim=1; end
    ax=nexttile(tl,j); imagesc(ax,T,[0 lim]); axis(ax,'image'); colorbar(ax);
    title(ax,['True ',defs{j,1}]);
    ax=nexttile(tl,4+j); imagesc(ax,E,[0 lim]); axis(ax,'image'); colorbar(ax);
    title(ax,['Estimated ',defs{j,1}]);
end
sgtitle(tl,sprintf('%s | %s | continuous interaction recovery',method,criterion), ...
    'Interpreter','none'); exportgraphics(f,[out,'_continuous_components.png'],'Resolution',180); close(f);

f=figure('Visible','off','Color','w','Position',[50 50 1500 720]);
tl=tiledlayout(f,2,4,'TileSpacing','compact','Padding','compact');
for j=1:4
    mask=t.truth_masks.(defs{j,3}); T=logical(t.truth_support_structural.(defs{j,2}));
    x=max(0,A.*mask); mx=max(x(mask)); if mx>0, x=x/mx; end
    E=(x>tau.(defs{j,4}))&mask;
    ax=nexttile(tl,j); imagesc(ax,T); axis(ax,'image'); title(ax,['True ',defs{j,1}]);
    ax=nexttile(tl,4+j); imagesc(ax,E); axis(ax,'image'); title(ax,['Estimated ',defs{j,1}]);
end
sgtitle(tl,sprintf('%s | %s | binary supports with frozen validation thresholds', ...
    method,criterion),'Interpreter','none');
exportgraphics(f,[out,'_binary_components.png'],'Resolution',180); close(f);
writematrix(A,[out,'_A_PRIMARY_CONTINUOUS.csv']);
end

function plot_mosaic_modes(e,t,out)
perm=align_modes(e.B,t.B_true); Bh=e.B(:,perm); Gh=e.Gamma(:,perm); Ah=e.A_state(:,:,perm);
writematrix(Bh,[out,'_B_ESTIMATED_ALIGNED.csv']);
writematrix(Gh,[out,'_GAMMA_ESTIMATED_ALIGNED.csv']);
for m=1:size(Ah,3), writematrix(Ah(:,:,m),sprintf('%s_A_MODE_%d_ESTIMATED.csv',out,m)); end
f=figure('Visible','off','Color','w','Position',[100 100 1250 650]);
tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); imagesc(ax,t.B_true); colorbar(ax); title(ax,'True B'); xlabel(ax,'mode'); ylabel(ax,'view');
ax=nexttile(tl); imagesc(ax,Bh); colorbar(ax); title(ax,'Estimated B (aligned)'); xlabel(ax,'mode'); ylabel(ax,'view');
ax=nexttile(tl,[1 2]); hold(ax,'on');
for m=1:size(Bh,2)
    plot(ax,1:t.K,t.B_true(:,m),'-o','LineWidth',1.8,'DisplayName',sprintf('true mode %d',m));
    plot(ax,1:t.K,Bh(:,m),'--s','LineWidth',1.5,'DisplayName',sprintf('estimated mode %d',m));
end
grid(ax,'on'); xlabel(ax,'view'); ylabel(ax,'simplex profile weight'); legend(ax,'Location','best');
sgtitle(tl,'MOSAIC view-profile recovery'); exportgraphics(f,[out,'_B_RECOVERY.png'],'Resolution',190); close(f);

r=size(Ah,3); f=figure('Visible','off','Color','w','Position',[50 50 450*r 760]);
tl=tiledlayout(f,2,r,'TileSpacing','compact','Padding','compact');
for m=1:r
    lim=max([t.A_state_true(:,:,m),Ah(:,:,m)],[],'all'); if lim<=0, lim=1; end
    ax=nexttile(tl,m); imagesc(ax,t.A_state_true(:,:,m),[0 lim]); axis(ax,'image'); colorbar(ax); title(ax,sprintf('True A_%d',m));
    ax=nexttile(tl,r+m); imagesc(ax,Ah(:,:,m),[0 lim]); axis(ax,'image'); colorbar(ax); title(ax,sprintf('Estimated A_%d',m));
end
sgtitle(tl,'MOSAIC physical modal patterns'); exportgraphics(f,[out,'_A_MODES.png'],'Resolution',190); close(f);

f=figure('Visible','off','Color','w','Position',[100 100 950 430]);
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); imagesc(ax,t.Gamma_true); colorbar(ax); title(ax,'True Gamma'); xlabel(ax,'mode'); ylabel(ax,'node');
ax=nexttile(tl); imagesc(ax,Gh); colorbar(ax); title(ax,'Estimated Gamma (aligned)'); xlabel(ax,'mode'); ylabel(ax,'node');
exportgraphics(f,[out,'_GAMMA_RECOVERY.png'],'Resolution',190); close(f);
end

function p=align_modes(B,T)
r=size(B,2); P=perms(1:r); loss=inf(size(P,1),1);
for j=1:size(P,1), loss(j)=norm(B(:,P(j,:))-T,'fro'); end
[~,q]=min(loss); p=P(q,:);
end

function plot_zw_primitive(e,t,out)
f=figure('Visible','off','Color','w','Position',[100 100 1250 420]);
tl=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
M={t.A_supra_true,e.A_primitive,e.A_effective}; names={'True supra','ZW primitive supra','ZW effective interaction'};
for j=1:3, ax=nexttile(tl); imagesc(ax,M{j}); axis(ax,'image'); colorbar(ax); title(ax,names{j}); end
exportgraphics(f,[out,'_ZW_PRIMITIVE_VS_EFFECTIVE.png'],'Resolution',180); close(f);
writematrix(e.A_primitive,[out,'_ZW_PRIMITIVE.csv']);
writematrix(e.A_effective,[out,'_ZW_EFFECTIVE.csv']);
end
