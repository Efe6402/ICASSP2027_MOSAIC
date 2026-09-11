function [truths,index,D]=prepare_density_bank(cfg)
%PREPARE_DENSITY_BANK Create exact-count nested truths and signal banks.
if ~isfolder(cfg.data_dir), mkdir(cfg.data_dir); end
manifest=fullfile(cfg.data_dir,'DENSITY_LEVELS.csv');
indexfile=fullfile(cfg.data_dir,'BANK_INDEX.csv');
bundle=fullfile(cfg.data_dir,'DENSITY_TRUTHS.mat');
if ~cfg.overwrite_data && isfile(manifest) && isfile(indexfile) && isfile(bundle)
    z=load(bundle,'truths'); truths=z.truths;
    D=readtable(manifest,'TextType','string');
    index=readtable(indexfile,'TextType','string');
    if valid_cached(truths,index,D,cfg), return; end
end

z=load(cfg.source_truth);
if isfield(z,'truth'), base=z.truth; elseif isfield(z,'data'), base=z.data;
else, error('Source truth MAT must contain truth or data.'); end
truths=cell(height_levels(cfg),1); rows=table();
truths{1}=normalize_truth(base,cfg,1);

source_counts=modal_edge_counts(truths{1}.A_state_true);
assert(isequal(source_counts,cfg.target_mode_edges(1,:)), ...
    ['The first row of target_mode_edges must equal the embedded source ', ...
     'counts. Source=%s, configured=%s.'],mat2str(source_counts), ...
     mat2str(cfg.target_mode_edges(1,:)));

old=rng; cleaner=onCleanup(@()rng(old));
rng(cfg.topology_seed+777,'twister');
priorities=cell(truths{1}.r,1);
for m=1:truths{1}.r
    P=rand(truths{1}.n); P=triu(P,1); priorities{m}=P;
end
truths{1}.density_priority=priorities;
for d=2:height_levels(cfg)
    truths{d}=expand_truth_exact(truths{d-1},cfg.target_mode_edges(d,:), ...
        cfg.q_in_levels(d,:),cfg.q_out_levels(d),priorities,cfg,d);
    assert(all(truths{d-1}.A_state_true(:)<=truths{d}.A_state_true(:)), ...
        'Physical supports are not nested at density level %d.',d);
end

for d=1:numel(truths)
    truths{d}.density_id=d;
    truths{d}.q_in=cfg.q_in_levels(d,:); truths{d}.q_out=cfg.q_out_levels(d);
    counts=modal_edge_counts(truths{d}.A_state_true);
    assert(isequal(counts,cfg.target_mode_edges(d,:)), ...
        'Exact modal-edge target was not attained at density %d.',d);
    [within_counts,outside_counts]=stratified_counts(truths{d});
    possible=nchoosek(truths{d}.n,2);
    target_density=mean(cfg.target_mode_edges(d,:))/possible;
    realized=mean(counts)/possible;
    rows=[rows;table(d,cfg.q_in_levels(d,1),cfg.q_in_levels(d,2), ...
        cfg.q_out_levels(d),cfg.contrast_odds_ratio,target_density,realized, ...
        cfg.target_mode_edges(d,1),cfg.target_mode_edges(d,2), ...
        counts(1),counts(2),within_counts(1),within_counts(2), ...
        outside_counts(1),outside_counts(2), ...
        'VariableNames',{'density_id','q_in_mode1','q_in_mode2','q_out', ...
        'contrast_odds_ratio','expected_mean_modal_density', ...
        'realized_mean_modal_density','target_mode1_edges','target_mode2_edges', ...
        'mode1_edges','mode2_edges','mode1_within_edges','mode2_within_edges', ...
        'mode1_background_edges','mode2_background_edges'})]; %#ok<AGROW>
    truth=truths{d};
    ddir=fullfile(cfg.data_dir,sprintf('density_%02d',d));
    if ~isfolder(ddir), mkdir(ddir); end
    save(fullfile(ddir,'TRUTH.mat'),'truth','-v7.3');
end
D=rows; writetable(D,manifest); save(bundle,'truths','-v7.3');
index=generate_signal_records(truths,cfg); writetable(index,indexfile);
end

function n=height_levels(cfg), n=size(cfg.target_mode_edges,1); end

function tf=valid_cached(truths,index,D,cfg)
nlev=height_levels(cfg); expected=nlev*(cfg.validation_count+cfg.test_count);
needed={'target_mode1_edges','target_mode2_edges','mode1_edges','mode2_edges'};
tf=numel(truths)==nlev && height(D)==nlev && height(index)==expected && ...
    all(index.sample_count==cfg.signal_count) && all(ismember(needed,D.Properties.VariableNames));
if ~tf, return; end
cached=[D.target_mode1_edges,D.target_mode2_edges];
tf=isequal(double(cached),double(cfg.target_mode_edges)) && ...
    isequal(double([D.mode1_edges,D.mode2_edges]),double(cfg.target_mode_edges));
end

function t=normalize_truth(d,cfg,id)
need={'n','K','r','A_state_true','B_true','Gamma_true','A_supra_true', ...
    'A_physical_true','A_copy_true','truth_masks','truth_support_structural', ...
    'Sigma','Sigma_clean','cfg','community_nodes'};
for j=1:numel(need), assert(isfield(d,need{j}),'Missing source field %s.',need{j}); end
t=struct(); for j=1:numel(need), t.(need{j})=d.(need{j}); end
optional={'Gamma_base_true','copy_scale_applied','copy_supports','Theta_true', ...
    'edge_i','edge_j','W_contribution_realized','topology_diagnostics'};
for j=1:numel(optional), if isfield(d,optional{j}), t.(optional{j})=d.(optional{j}); end, end
t.N=t.n*t.K; t.density_id=id; t.source='embedded MOSAIC truth';
t.A_within_true=t.A_supra_true.*t.truth_masks.within;
t.A_cross_true=t.A_supra_true.*t.truth_masks.crossnode_crossview;
t.A_copy_masked_true=t.A_supra_true.*t.truth_masks.copy;
t.signal_count=cfg.signal_count;
end

function t=expand_truth_exact(prev,target_edges,qin,qout,priorities,cfg,id)
% Add an exact number of edges while approaching the fixed-odds allocation.
n=prev.n; A=prev.A_state_true>0; upper=triu(true(n),1);
for m=1:prev.r
    C=prev.community_nodes{m}; inside=false(n); inside(C,C)=true;
    inside=inside&upper; outside=upper&~inside; old=A(:,:,m)&upper;
    current=nnz(old); target=target_edges(m);
    assert(target>=current && target<=nnz(upper), ...
        'Invalid exact edge target for mode %d at density %d.',m,id);
    need=target-current;
    if need==0, continue; end

    mass_in=nnz(inside)*qin(m); mass_out=nnz(outside)*qout;
    desired_in=round(target*mass_in/max(mass_in+mass_out,eps));
    desired_in=min(max(desired_in,0),nnz(inside));
    current_in=nnz(old&inside);
    add_in=min(max(desired_in-current_in,0),need);
    add_out=need-add_in;
    avail_in=nnz(inside&~old); avail_out=nnz(outside&~old);
    if add_in>avail_in
        add_out=add_out+(add_in-avail_in); add_in=avail_in;
    end
    if add_out>avail_out
        add_in=add_in+(add_out-avail_out); add_out=avail_out;
    end
    assert(add_in<=avail_in && add_out<=avail_out && add_in+add_out==need, ...
        'Unable to allocate the requested exact edge count.');

    U=old;
    U=activate_lowest_priority(U,inside&~old,priorities{m},add_in);
    U=activate_lowest_priority(U,outside&~old,priorities{m},add_out);
    A(:,:,m)=U|U';
end

t=prev; t.A_state_true=double(A); t.density_id=id;
t.density_priority=priorities;
[ei,ej]=find(upper); t.edge_i=ei; t.edge_j=ej;
t.Theta_true=zeros(numel(ei),prev.r);
for m=1:prev.r
    Am=t.A_state_true(:,:,m);
    t.Theta_true(:,m)=Am(sub2ind([n,n],ei,ej));
end

% The numerical copy coefficients are deliberately held fixed. Only the
% physical modal supports change along this density path.
t.Gamma_true=prev.Gamma_true;
if isfield(prev,'Gamma_base_true'), t.Gamma_base_true=prev.Gamma_base_true; end
if isfield(prev,'copy_scale_applied'), t.copy_scale_applied=prev.copy_scale_applied; end
[t.A_physical_true,t.A_copy_true]=supra_parts( ...
    t.A_state_true,t.Gamma_true,t.B_true);
t.A_supra_true=zero_sym(t.A_physical_true+t.A_copy_true);
t.truth_support_structural=struct( ...
    'all',t.A_supra_true>0, ...
    'within',(t.A_supra_true>0)&t.truth_masks.within, ...
    'crossnode_crossview',(t.A_supra_true>0)&t.truth_masks.crossnode_crossview, ...
    'copy',(t.A_supra_true>0)&t.truth_masks.copy);
t.A_within_true=t.A_supra_true.*t.truth_masks.within;
t.A_cross_true=t.A_supra_true.*t.truth_masks.crossnode_crossview;
t.A_copy_masked_true=t.A_supra_true.*t.truth_masks.copy;
t.cfg.q_in=qin; t.cfg.q_out=qout; t.cfg.p=cfg.signal_count;
[t.Sigma_clean,t.Sigma]=signal_covariance(t.A_supra_true,t.cfg);
t.source='exact-count nested density expansion from embedded MOSAIC truth';
end

function U=activate_lowest_priority(U,available,P,count)
if count==0, return; end
idx=find(available); [~,order]=sort(P(idx),'ascend');
assert(count<=numel(order),'Not enough candidate edges to activate.');
U(idx(order(1:count)))=true;
end

function counts=modal_edge_counts(A)
r=size(A,3); counts=zeros(1,r);
for m=1:r, counts(m)=nnz(triu(A(:,:,m)>0,1)); end
end

function [within_counts,outside_counts]=stratified_counts(t)
upper=triu(true(t.n),1); within_counts=zeros(1,t.r); outside_counts=zeros(1,t.r);
for m=1:t.r
    C=t.community_nodes{m}; inside=false(t.n); inside(C,C)=true; inside=inside&upper;
    A=t.A_state_true(:,:,m)>0;
    within_counts(m)=nnz(A&inside); outside_counts(m)=nnz(A&(upper&~inside));
end
end

function [Aphys,Acopy]=supra_parts(Astate,Gamma,B)
[n,~,r]=size(Astate); K=size(B,1); N=n*K; Aphys=zeros(N); Acopy=zeros(N);
for m=1:r
    bb=B(:,m)*B(:,m)';
    Aphys=Aphys+kron(bb,Astate(:,:,m));
    Acopy=Acopy+kron(bb-diag(diag(bb)),diag(Gamma(:,m)));
end
Aphys=zero_sym(Aphys); Acopy=zero_sym(Acopy);
end

function A=zero_sym(A), A=.5*(A+A'); A(1:size(A,1)+1:end)=0; end

function [Sclean,Sigma]=signal_covariance(A,c)
deg=sum(A,2); L=diag(deg)-A; lm=max(real(eig(L))); if lm>0, L=L/lm; end
[V,D]=eig(.5*(L+L')); ell=max(real(diag(D)),0);
if strcmpi(c.signal_model,'heat_diffusion'), h=exp(-c.heat_time*ell);
else, h=(1+c.filter_strength*ell).^(-c.filter_order); end
H=V*diag(h)*V'; Sclean=c.signal_variance*(H*H'); Sclean=.5*(Sclean+Sclean');
Sigma=Sclean+c.observation_noise_variance*eye(size(A,1)); Sigma=.5*(Sigma+Sigma');
end

function index=generate_signal_records(truths,cfg)
split_col=strings(0,1); density_col=[]; trial_col=[]; seed_col=[]; file_col=strings(0,1);
splits=["validation","test"]; counts=[cfg.validation_count,cfg.test_count];
for d=1:numel(truths)
    t=truths{d}; [V,D]=eig(.5*(t.Sigma+t.Sigma'));
    F=real(V*diag(sqrt(max(real(diag(D)),0))));
    ddir=fullfile(cfg.data_dir,sprintf('density_%02d',d));
    for si=1:2
        for trial=1:counts(si)
            seed=cfg.seed_base+d*1000000+si*100000+trial;
            old=rng; clean=onCleanup(@()rng(old)); rng(seed,'twister');
            Y=F*randn(t.N,cfg.signal_count); record=make_record(Y,t,cfg, ...
                char(splits(si)),trial,seed);
            file=fullfile(ddir,sprintf('%s_trial%02d.mat',splits(si),trial));
            save(file,'record','-v7.3'); clear clean
            split_col(end+1,1)=splits(si); density_col(end+1,1)=d; %#ok<AGROW>
            trial_col(end+1,1)=trial; seed_col(end+1,1)=seed; file_col(end+1,1)=file; %#ok<AGROW>
        end
    end
end
index=table(split_col,density_col,repmat(cfg.signal_count,numel(density_col),1), ...
    trial_col,seed_col,file_col,'VariableNames', ...
    {'split','density_id','sample_count','trial','seed','file'});
end

function record=make_record(Y,t,cfg,split,trial,seed)
p=cfg.signal_count; Cy=(Y*Y')/p; d=diag(Cy); S=max(d+d'-2*Cy,0); S=zero_sym(S);
Sp=zeros(t.n); Sq=zeros(t.K);
for s=1:p
    Xi=reshape(Y(:,s),t.n,t.K); Sp=Sp+cov(Xi'); Sq=Sq+cov(Xi);
end
record=struct('split',split,'density_id',t.density_id,'sample_count',p, ...
    'trial',trial,'seed',seed,'n',t.n,'K',t.K,'S',S, ...
    'Cy_uncentered',Cy,'Sp',Sp/p,'Sq',Sq/p);
end
