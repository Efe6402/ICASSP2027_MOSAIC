%% Combined rank-one supra-F1 figure for the three experiments
clearvars; close all; clc;
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'conference_figures'); if ~isfolder(out), mkdir(out); end

pfile=latest_result(fullfile(root,'01_sample_size_and_model_selection', ...
    'results','p5000_fixed_transfer'),'ALL_TEST_SUMMARY.csv');
dfile=latest_result(fullfile(root,'02_density_sweep','results','density_sweep'), ...
    'TEST_SUMMARY.csv');
rfile=latest_result(fullfile(root,'03_view_profile_rho_sweep','results', ...
    'rho_sweep_p05000'),'TEST_SUMMARY.csv');

Tp=readtable(pfile,'TextType','string');
Tp=Tp(Tp.selection_criterion=="f1_component"&Tp.selection_rank==1& ...
    Tp.threshold_policy=="anchor_frozen",:);
Td=readtable(dfile,'TextType','string');
Td=Td(Td.selection_criterion=="f1_component"&Td.selection_rank==1& ...
    Td.threshold_policy=="p5000_validation_fixed",:);
Tr=readtable(rfile,'TextType','string');
Tr=Tr(Tr.selection_criterion=="f1_component"&Tr.selection_rank==1& ...
    Tr.threshold_policy=="p5000_validation_fixed",:);
assert_grid(Tp,'evaluation_sample_count','sample-size');
assert_grid(Td,'realized_density','density'); assert_grid(Tr,'rho','rho');

methods=["MOSAIC","PGL2021","ZW2024"];
labels=["MOSAIC","PGL2021","ZW2024"];
colors=[0 .447 .741;.85 .325 .098;.929 .694 .125];
markers={'o','s','^'};
fig=figure('Color','w','Position',[80 60 1650 1080]);
width=.385; top_y=.60; top_h=.335; gap=.145; bottom_h=.30;
ax1=axes(fig,'Position',[.065 top_y width top_h]);
h=panel(ax1,Tp,'evaluation_sample_count','test_f1_supra_mean', ...
    'test_f1_supra_std','Number of graph signals $p$', ...
    '(i) Supra-adjacency F1 vs. $p$',methods,colors,markers);
set(ax1,'XScale','log'); xticks(ax1,[20 100 1000 10000]);
ax2=axes(fig,'Position',[.555 top_y width top_h]);
panel(ax2,Td,'realized_density','test_f1_supra_mean','test_f1_supra_std', ...
    'Realized mean modal edge density','(ii) Supra-adjacency F1 vs. density', ...
    methods,colors,markers);
bottom_y=top_y-gap-bottom_h;
ax3=axes(fig,'Position',[(1-width)/2 bottom_y width bottom_h]);
panel(ax3,Tr,'rho','f1_supra_mean','f1_supra_std','$\rho$', ...
    '(iii) Supra-adjacency F1 vs. $\rho$',methods,colors,markers);
lgd=legend(ax3,h,cellstr(labels),'Orientation','horizontal','FontSize',16, ...
    'Box','off','Interpreter','none','Location','none');
lgd.Units='normalized'; lgd.Position=[.325 .02 .35 .05];
drawnow;
name='TOP1_FIXED_VALIDATION_TAU_SUPRA_F1_THREE_SWEEPS';
exportgraphics(fig,fullfile(out,[name,'.png']),'Resolution',400);
exportgraphics(fig,fullfile(out,[name,'.pdf']),'ContentType','vector');
savefig(fig,fullfile(out,[name,'.fig']));
fprintf('\nSaved:\n%s\n',fullfile(out,[name,'.png']));

function file=latest_result(base,name)
assert(isfolder(base),'Results directory not found: %s',base);
D=dir(fullfile(base,'run_*')); D=D([D.isdir]);
[~,order]=sort([D.datenum],'descend'); D=D(order);
file='';
for k=1:numel(D)
    candidate=fullfile(D(k).folder,D(k).name,name);
    complete=fullfile(D(k).folder,D(k).name,'COMPLETE_RESULTS.mat');
    if isfile(candidate)&&isfile(complete), file=candidate; break; end
end
assert(~isempty(file),'No completed result containing %s was found under %s.',name,base);
end

function assert_grid(T,xname,label)
methods=["MOSAIC","PGL2021","ZW2024"]; ref=[];
for k=1:numel(methods)
    Q=T(upper(T.method)==methods(k),:); assert(~isempty(Q),'Missing %s %s rows.',methods(k),label);
    x=sort(double(Q.(xname))); assert(numel(x)==numel(unique(x)));
    if isempty(ref), ref=x; else, assert(isequal(ref,x),'Inconsistent %s grid.',label); end
end
end

function h=panel(ax,T,xname,mn,sd,xlabel_text,title_text,methods,colors,markers)
hold(ax,'on'); grid(ax,'on'); box(ax,'on'); h=gobjects(numel(methods),1);
for k=1:numel(methods)
    Q=T(upper(string(T.method))==methods(k),:); Q=sortrows(Q,xname);
    h(k)=errorbar(ax,double(Q.(xname)),double(Q.(mn)),double(Q.(sd)),'-', ...
        'Color',colors(k,:),'Marker',markers{k},'MarkerSize',8, ...
        'MarkerFaceColor','w','LineWidth',2.3,'CapSize',8);
end
xlabel(ax,xlabel_text,'Interpreter','latex','FontSize',25);
ylabel(ax,'Supra F1','Interpreter','latex','FontSize',25);
title(ax,title_text,'Interpreter','latex','FontWeight','bold','FontSize',20);
ylim(ax,[0 1.02]); set(ax,'FontName','Helvetica','FontSize',17, ...
    'LineWidth',1.15,'TickDir','out','GridAlpha',.18);
end
