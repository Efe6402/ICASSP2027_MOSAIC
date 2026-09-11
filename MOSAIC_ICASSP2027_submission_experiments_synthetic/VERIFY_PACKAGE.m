%% Static integrity checks for the submission package
clearvars; clc;
root=fileparts(mfilename('fullpath'));
packages=["01_sample_size_and_model_selection", ...
    "02_density_sweep","03_view_profile_rho_sweep"];
methods=["MOSAIC","PGL2021","ZW2024"];
expected_ids=struct('MOSAIC',[407 420 421], ...
    'PGL2021',[269 294 319],'ZW2024',[593 482 592]);

for p=packages
    pkg=fullfile(root,p); assert(isfolder(pkg),'Missing subpackage %s.',p);
    truth=fullfile(pkg,'data','source_truth','MOSAIC_PLANTED_BLOCK_TRUTH_ONLY.mat');
    assert(isfile(truth),'Missing embedded truth in %s.',p);
    for method=methods
        file=fullfile(pkg,'selected_configurations','p_05000', ...
            method+"_F1_COMPONENT_TOP3.csv");
        T=readtable(file,'TextType','string');
        assert(height(T)==3&&all(T.selection_criterion=="f1_component"));
        assert(isequal(double(T.candidate_id(:))',expected_ids.(char(method))));
        assert(isequal(double(T.selection_rank(:))',[1 2 3]));
        assert(all(isfinite(double(T{:,{'tau_supra','tau_within','tau_cross','tau_copy'}})),'all'));
    end
end

solver_rel={ ...
    fullfile('solver','mosaic','solve_mosaic_crossview.m'), ...
    fullfile('solver','pgl2021','Learn_PGL.m'), ...
    fullfile('solver','pgl2021','PGL_solver.m'), ...
    fullfile('solver','zw2024','solve_zhang_wai_2024_paper.m')};
for j=1:numel(solver_rel)
    reference=fileread(fullfile(root,packages(1),solver_rel{j}));
    for p=packages(2:end)
        assert(strcmp(reference,fileread(fullfile(root,p,solver_rel{j}))), ...
            'Solver contents differ across subpackages: %s',solver_rel{j});
    end
end

sources=dir(fullfile(root,'**','*.m'));
for k=1:numel(sources)
    text=fileread(fullfile(sources(k).folder,sources(k).name));
    marker=['/','Users','/'];
    assert(~contains(text,marker),'Absolute user path found in %s.',sources(k).name);
end
fprintf('Package verification passed: %d subpackages, %d selected configurations.\n', ...
    numel(packages),numel(packages)*numel(methods)*3);
