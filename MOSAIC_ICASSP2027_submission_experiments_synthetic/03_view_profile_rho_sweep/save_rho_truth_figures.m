function save_rho_truth_figures(cfg)
%SAVE_RHO_TRUTH_FIGURES Export truth matrices and paper-style diagnostics.
out=fullfile(cfg.output_dir,'TRUTH_BY_RHO'); if ~isfolder(out), mkdir(out); end
base=load_default_truth(cfg); T=cell(numel(cfg.rho_values),1); cmax=0;
for j=1:numel(T)
    T{j}=build_rho_truth(base,cfg.rho_values(j));
    cmax=max(cmax,max(T{j}.A_supra_true,[],'all'));
end
profile_rows=table();
for j=1:numel(T)
    t=T{j}; folder=fullfile(out,sprintf('rho_%02d',j));
    if ~isfolder(folder), mkdir(folder); end
    truth=t; save(fullfile(folder,'TRUTH.mat'),'truth','-v7.3');
    writematrix(t.W_true,fullfile(folder,'W_VIEW_CONTRIBUTIONS.csv'));
    writematrix(t.B_true,fullfile(folder,'B_COLUMN_SIMPLEX.csv'));
    writematrix(t.A_supra_true,fullfile(folder,'A_SUPRA_CONTINUOUS.csv'));
    writematrix(double(t.truth_support_structural.all),fullfile(folder,'SUPPORT_SUPRA.csv'));
    writematrix(double(t.truth_support_structural.within),fullfile(folder,'SUPPORT_WITHIN.csv'));
    writematrix(double(t.truth_support_structural.crossnode_crossview),fullfile(folder,'SUPPORT_CROSS.csv'));
    writematrix(double(t.truth_support_structural.copy),fullfile(folder,'SUPPORT_COPY.csv'));
    for k=1:t.K
        profile_rows=[profile_rows;table(j,t.rho,k,t.W_true(k,1),t.W_true(k,2), ...
            t.B_true(k,1),t.B_true(k,2),'VariableNames',{'level','rho','view', ...
            'W_mode1','W_mode2','B_mode1','B_mode2'})]; %#ok<AGROW>
    end
    save_truth_panel(t,folder,cmax,j);
    save_profile_panel(t,folder,j);
end
writetable(profile_rows,fullfile(out,'ALL_RHO_VIEW_PROFILES.csv'));
end

function save_truth_panel(t,folder,cmax,lev)
f=figure('Visible','off','Color','w','Position',[50 50 1500 920]);
tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('Generated truth at \\rho=%.5g (level %d)',t.rho,lev), ...
    'FontWeight','bold');
draw(nexttile,t.A_supra_true,sprintf('Continuous weighted supra-adjacency (%d x %d)',t.N,t.N),[0 cmax],t.n);
draw(nexttile,double(t.truth_support_structural.within),'Within-view structural support',[0 1],t.n);
draw(nexttile,double(t.truth_support_structural.crossnode_crossview),'Cross-node/cross-view support',[0 1],t.n);
draw(nexttile,double(t.truth_support_structural.copy),'Same-node copy support',[0 1],t.n);
exportgraphics(f,fullfile(folder,'TRUTH_SUPRA_AND_COMPONENT_SUPPORTS.png'),'Resolution',220);
exportgraphics(f,fullfile(folder,'TRUTH_SUPRA_AND_COMPONENT_SUPPORTS.pdf'),'ContentType','vector');
close(f);
end

function save_profile_panel(t,folder,lev)
f=figure('Visible','off','Color','w','Position',[100 100 1100 430]);
tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('View-profile construction at \\rho=%.5g (level %d)',t.rho,lev));
ax=nexttile; imagesc(ax,t.W_true); axis(ax,'image'); colorbar(ax); clim(ax,[0 1]);
title(ax,'Contribution fractions W'); xlabel(ax,'mode'); ylabel(ax,'view');
xticks(ax,1:2); yticks(ax,1:4);
ax=nexttile; imagesc(ax,t.B_true); axis(ax,'image'); colorbar(ax); clim(ax,[0 max(t.B_true,[],'all')]);
title(ax,'Column-simplex profiles B'); xlabel(ax,'mode'); ylabel(ax,'view');
xticks(ax,1:2); yticks(ax,1:4);
exportgraphics(f,fullfile(folder,'W_AND_B_PROFILES.png'),'Resolution',220);
exportgraphics(f,fullfile(folder,'W_AND_B_PROFILES.pdf'),'ContentType','vector');
close(f);
end

function draw(ax,A,label,lims,n)
imagesc(ax,A); axis(ax,'image'); title(ax,label,'FontWeight','bold');
xlabel(ax,'supra-node'); ylabel(ax,'supra-node'); colorbar(ax); clim(ax,lims);
hold(ax,'on'); N=size(A,1);
for q=n+.5:n:N
    xline(ax,q,'w-','LineWidth',.6); yline(ax,q,'w-','LineWidth',.6);
end
hold(ax,'off');
end
