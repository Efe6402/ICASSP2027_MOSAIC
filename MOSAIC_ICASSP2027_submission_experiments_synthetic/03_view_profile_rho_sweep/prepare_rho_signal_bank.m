function [truths,index,levels]=prepare_rho_signal_bank(cfg)
%PREPARE_RHO_SIGNAL_BANK Paired Gaussian records for every rho.
if ~isfolder(cfg.data_dir), mkdir(cfg.data_dir); end
indexfile=fullfile(cfg.data_dir,'BANK_INDEX.csv');
levelsfile=fullfile(cfg.data_dir,'RHO_LEVELS.csv');
truthfile=fullfile(cfg.data_dir,'RHO_TRUTHS.mat');
expected=numel(cfg.rho_values)*(cfg.validation_count+cfg.test_count);
if ~cfg.overwrite_data&&isfile(indexfile)&&isfile(levelsfile)&&isfile(truthfile)
    index=readtable(indexfile,'TextType','string'); levels=readtable(levelsfile);
    z=load(truthfile,'truths'); truths=z.truths;
    good=height(index)==expected&&numel(truths)==numel(cfg.rho_values)&& ...
        all(index.sample_count==cfg.anchor_p)&& ...
        isequal(double(levels.rho(:)'),double(cfg.rho_values(:)'));
    if good, return; end
end

base=load_default_truth(cfg); truths=cell(numel(cfg.rho_values),1);
levels=table((1:numel(cfg.rho_values))',cfg.rho_values(:), ...
    abs(cfg.rho_values(:)-cfg.rho_reference)<1e-12, ...
    'VariableNames',{'level','rho','is_source_reference'});
for lev=1:numel(cfg.rho_values)
    truths{lev}=build_rho_truth(base,cfg.rho_values(lev));
end

splits=["validation","test"];
counts=[cfg.validation_count,cfg.test_count]; rows=table();
for si=1:2
    for trial=1:counts(si)
        % The innovation seeds pair realizations across rho values.
        seed=cfg.seed_base+si*100000+trial;
        old=rng; clean=onCleanup(@()rng(old)); rng(seed,'twister');
        Z=randn(base.N,cfg.anchor_p);
        for lev=1:numel(truths)
            truth=truths{lev};
            [V,D]=eig(.5*(truth.Sigma+truth.Sigma'));
            F=real(V*diag(sqrt(max(real(diag(D)),0))));
            Y=F*Z;
            meta=struct('split',char(splits(si)),'level',lev,'trial',trial, ...
                'seed',seed,'innovation_seed',seed,'rho',truth.rho, ...
                'source_anchor_p',cfg.anchor_p);
            record=make_signal_record(Y,truth,meta);
            folder=fullfile(cfg.data_dir,sprintf('rho_%02d',lev));
            if ~isfolder(folder), mkdir(folder); end
            file=fullfile(folder,sprintf('%s_trial%02d.mat',splits(si),trial));
            save(file,'record','truth','-v7.3');
            rows=[rows;table(splits(si),lev,truth.rho,cfg.anchor_p,trial,seed, ...
                string(file),'VariableNames',{'split','level','rho', ...
                'sample_count','trial','innovation_seed','file'})]; %#ok<AGROW>
        end
        clear clean
    end
end
index=rows; writetable(index,indexfile); writetable(levels,levelsfile);
save(truthfile,'truths','-v7.3');
end
