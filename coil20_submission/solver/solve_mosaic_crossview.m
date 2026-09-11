function result = solve_mosaic_crossview(observed, opts)
%SOLVE_MOSAIC_CROSSVIEW Three-block PALM solver for no-mass MOSAIC.
%
% Variables
%   Theta >= 0 : E-by-r cross-node modal edge coefficients
%   Gamma >= 0 : n-by-r same-node copy coefficients
%   B           : K-by-r, each column in the probability simplex
%
% The objective is the simplified no-prescribed-mass MOSAIC objective:
% cross-node and copy smoothness, l1/l2 coefficient penalties, two
% complementary smoothed log-degree penalties, and modal-profile
% diversity.  The two penalized degrees are
%   D_intra = (H*Theta)*(B.^2)'
%   D_inter = (H*Theta+Gamma)*(B.*(1-B))'.
% D_supra=D_intra+D_inter and d_view=sum(D_inter,1)' are derived
% diagnostics and are not penalized again.
%
% Important: epsilon-smoothed logarithms are penalties, not strict barriers.
% With beta_a,beta_c>0, positive epsilons, positive line-search margins,
% sufficient-decrease backtracking, and bounded simplex B, PALM theory
% supports convergence of the iterate sequence to a critical point (not a
% global optimum). Initial profiles are random interior-simplex points and
% contain no planted modal-center information.

if nargin < 2, opts = struct(); end
opts = defaults(opts);
[problem, observed] = prepare_problem(observed, opts);
if strcmpi(opts.operation, 'gradient_check')
    result = gradient_check(problem, opts);
    return
end

old_rng = rng;
cleanup = onCleanup(@() rng(old_rng));
rng(opts.seed, 'twister');

records = cell(opts.n_restarts,1);
bestJ = inf;
best = 0;
for s = 1:opts.n_restarts
    [Theta,Gamma,B] = initialize(problem,opts,s);
    records{s} = one_restart(Theta,Gamma,B,problem,opts,s);
    if isfinite(records{s}.objective) && records{s}.objective < bestJ
        bestJ = records{s}.objective;
        best = s;
    end
    if opts.verbose
        fprintf('MOSAIC restart %d/%d: J=% .6e, residual=%.3e, %s\n', ...
            s,opts.n_restarts,records{s}.objective, ...
            records{s}.stationarity_combined,records{s}.status);
    end
end
if best == 0, error('Every MOSAIC restart failed.'); end
result = records{best};
result.best_restart = best;
result.all_restarts = records;
result.options = opts;
result.cost_scale = problem.cost_scale;
result.observed = observed;
end

function rec = one_restart(Theta,Gamma,B,P,o,restart)
Lth=o.L_theta; Lga=o.L_gamma; LB=o.L_B;
B_initial=B;
hist=nan(o.max_iter+1,1);
res_hist=nan(o.max_iter+1,3);
[J,Gth,~,~,~]=objective_gradient(Theta,Gamma,B,P,o);
hist(1)=J;
status='maximum_iterations';

for it=1:o.max_iter
    J0=J;

    % Theta update.
    [Theta_new,Lth,ok]=backtrack_nonnegative(Theta,Gth,Lth, ...
        @(X)objective_only(X,Gamma,B,P,o),J,o,o.sigma_theta);
    if ~ok, status='theta_line_search_failure'; break; end
    Theta=Theta_new;

    % Gamma update at the new Theta.
    [Jmid,~,Gga,~,~]=objective_gradient(Theta,Gamma,B,P,o);
    [Gamma_new,Lga,ok]=backtrack_nonnegative(Gamma,Gga,Lga, ...
        @(X)objective_only(Theta,X,B,P,o),Jmid,o,o.sigma_gamma);
    if ~ok, status='gamma_line_search_failure'; break; end
    Gamma=Gamma_new;

    % B update at the new Theta and Gamma.
    [Jmid,~,~,GB,~]=objective_gradient(Theta,Gamma,B,P,o);
    [B_new,LB,ok]=backtrack_simplex(B,GB,LB, ...
        @(X)objective_only(Theta,Gamma,X,P,o),Jmid,o,o.sigma_B);
    if ~ok, status='B_line_search_failure'; break; end
    B=B_new;

    [J,Gth,Gga,GB,~]=objective_gradient(Theta,Gamma,B,P,o);
    hist(it+1)=J;
    rth=projected_residual_nonnegative(Theta,Gth,Lth);
    rga=projected_residual_nonnegative(Gamma,Gga,Lga);
    rB=projected_residual_simplex(B,GB,LB);
    res_hist(it+1,:)=[rth,rga,rB];
    relJ=abs(J-J0)/max(1,abs(J0));
    if max([rth,rga,rB]) <= o.stationarity_tol && relJ <= o.objective_tol
        status='converged';
        break
    end
    if ~isfinite(J), status='nonfinite_objective'; break; end
end

last=find(isfinite(hist),1,'last');
if isempty(last), last=1; end
hist=hist(1:last);
res_hist=res_hist(1:last,:);
[J,Gth,Gga,GB,aux]=objective_gradient(Theta,Gamma,B,P,o);
rth=projected_residual_nonnegative(Theta,Gth,Lth);
rga=projected_residual_nonnegative(Gamma,Gga,Lga);
rB=projected_residual_simplex(B,GB,LB);
graph=reconstruct(Theta,Gamma,B,P);

rec=graph;
rec.Theta=Theta; rec.Gamma=Gamma; rec.B=B;
rec.initial_B=B_initial;
rec.objective=J; rec.objective_components=aux.components;
rec.status=status; rec.restart=restart; rec.iterations=last-1;
rec.stationarity_theta=rth;
rec.stationarity_gamma=rga;
rec.stationarity_B=rB;
rec.stationarity_combined=sqrt(rth^2+rga^2+rB^2);
rec.objective_history=hist;
rec.residual_history=res_hist;
rec.minimum_intra_degree=min(aux.D_intra(:));
rec.minimum_inter_degree=min(aux.D_inter(:));
rec.minimum_copy_degree=min(aux.D_copy(:));
rec.minimum_supra_degree=min(aux.D_supra(:));
rec.minimum_view_degree=min(aux.d_view(:));
rec.feasibility=struct( ...
    'theta_nonnegativity',max(0,-min(Theta(:))), ...
    'gamma_nonnegativity',max(0,-min(Gamma(:))), ...
    'B_nonnegativity',max(0,-min(B(:))), ...
    'B_simplex',max(abs(sum(B,1)-1)));
end

function [Xnew,L,ok]=backtrack_nonnegative(X,G,L,fun,J,o,sigma)
L=max([o.L_min,L/o.bt_relax,sigma*(1+1e-12)]);
ok=false;
for q=1:o.max_backtracks
    Xnew=max(X-G/L,0);
    d=Xnew-X;
    Jnew=fun(Xnew);
    upper=J+sum(G(:).*d(:))+0.5*(L-sigma)*sum(d(:).^2);
    if isfinite(Jnew) && Jnew <= upper+o.bt_tol
        ok=true; return
    end
    L=L*o.bt_factor;
end
Xnew=X;
end

function [Bnew,L,ok]=backtrack_simplex(B,G,L,fun,J,o,sigma)
L=max([o.L_min,L/o.bt_relax,sigma*(1+1e-12)]);
ok=false;
for q=1:o.max_backtracks
    Bnew=project_simplex_columns(B-G/L);
    d=Bnew-B;
    Jnew=fun(Bnew);
    upper=J+sum(G(:).*d(:))+0.5*(L-sigma)*sum(d(:).^2);
    if isfinite(Jnew) && Jnew <= upper+o.bt_tol
        ok=true; return
    end
    L=L*o.bt_factor;
end
Bnew=B;
end

function [J,Gth,Gga,GB,A]=objective_gradient(Theta,Gamma,B,P,o)
Q=P.H*Theta;
U=B.^2;
W=B.*(1-B);
sa=sum(Theta,1)';
sc=sum(Gamma,1)';
h=2*sa+sc;
D_intra=Q*U';
D_cross=Q*W';
D_copy=Gamma*W';
D_inter=D_cross+D_copy;
D_supra=D_intra+D_inter;
d_view=W*h;

Da=zeros(P.E,P.r); Dc=zeros(P.n,P.r);
Ta=zeros(P.K,P.K,P.r); Tc=zeros(P.K,P.K,P.r);
for m=1:P.r
    b=B(:,m);
    for e=1:P.E
        Se=P.S_edge(:,:,e);
        Da(e,m)=b'*Se*b;
        Ta(:,:,m)=Ta(:,:,m)+Theta(e,m)*Se;
    end
    for i=1:P.n
        Ri=P.R_copy(:,:,i);
        Dc(i,m)=0.5*o.eta_c*(b'*Ri*b);
        Tc(:,:,m)=Tc(:,:,m)+Gamma(i,m)*Ri;
    end
end

div=0;
for m=1:P.r
    for q=m+1:P.r, div=div+(B(:,m)'*B(:,q))^2; end
end
c_data=sum(Theta(:).*Da(:))+sum(Gamma(:).*Dc(:));
c_l1=o.alpha_a*sum(Theta(:))+o.alpha_c*sum(Gamma(:));
c_l2=0.5*o.beta_a*sum(Theta(:).^2)+0.5*o.beta_c*sum(Gamma(:).^2);
c_intra=-o.tau_intra*sum(log(D_intra(:)+o.epsilon_intra));
c_inter=-o.tau_inter*sum(log(D_inter(:)+o.epsilon_inter));
c_B=o.lambda_B*div;
J=c_data+c_l1+c_l2+c_intra+c_inter+c_B;

Uintra=1./(D_intra+o.epsilon_intra);
Uinter=1./(D_inter+o.epsilon_inter);
Gth=Da+o.alpha_a+o.beta_a*Theta ...
    -o.tau_intra*P.H'*(Uintra*U) ...
    -o.tau_inter*P.H'*(Uinter*W);
Gga=Dc+o.alpha_c+o.beta_c*Gamma ...
    -o.tau_inter*Uinter*W;

GB=zeros(P.K,P.r);
for m=1:P.r
    b=B(:,m);
    q_intra=sum((Q(:,m).*ones(1,P.K)).*Uintra,1)';
    q_inter=sum(((Q(:,m)+Gamma(:,m)).*ones(1,P.K)).*Uinter,1)';
    GB(:,m)=2*Ta(:,:,m)*b+o.eta_c*Tc(:,:,m)*b ...
        -2*o.tau_intra*b.*q_intra ...
        -o.tau_inter*(1-2*b).*q_inter;
    for q=1:P.r
        if q~=m
            GB(:,m)=GB(:,m)+2*o.lambda_B*(b'*B(:,q))*B(:,q);
        end
    end
end

A=struct();
A.D_intra=D_intra; A.D_cross=D_cross;
A.D_copy=D_copy; A.D_inter=D_inter;
A.D_supra=D_supra; A.d_view=d_view;
A.components=struct('data',c_data,'l1',c_l1,'l2',c_l2, ...
    'log_intra',c_intra,'log_inter',c_inter,'diversity',c_B);
end

function J=objective_only(Theta,Gamma,B,P,o)
[J,~,~,~,~]=objective_gradient(Theta,Gamma,B,P,o);
end

function [P,observed]=prepare_problem(observed,o)
required={'n','K','r'};
for q=1:numel(required)
    if ~isfield(observed,required{q}), error('observed.%s is required.',required{q}); end
end
n=observed.n; K=observed.K; r=observed.r;
[ii,jj]=find(triu(true(n),1)); E=numel(ii);
H=sparse([ii;jj],[(1:E)';(1:E)'],1,n,E);

if isfield(observed,'S_edge') && isfield(observed,'R_copy')
    Sedge=observed.S_edge; Rcopy=observed.R_copy;
elseif isfield(observed,'S') && isequal(size(observed.S),[n*K,n*K])
    Sedge=zeros(K,K,E); Rcopy=zeros(K,K,n);
    Sfull=0.5*(observed.S+observed.S');
    for e=1:E
        for k=1:K
            for ell=1:K
                Sedge(k,ell,e)=Sfull((k-1)*n+ii(e),(ell-1)*n+jj(e));
            end
        end
        Sedge(:,:,e)=0.5*(Sedge(:,:,e)+Sedge(:,:,e)');
    end
    for i=1:n
        for k=1:K
            for ell=1:K
                if k~=ell
                    Rcopy(k,ell,i)=Sfull((k-1)*n+i,(ell-1)*n+i);
                end
            end
        end
        Ri=0.5*(Rcopy(:,:,i)+Rcopy(:,:,i)');
        Ri(1:K+1:end)=0;
        Rcopy(:,:,i)=Ri;
    end
else
    if ~isfield(observed,'Y'), error('Provide observed.S or observed.Y.'); end
    Y=observed.Y;
    if size(Y,1)~=n*K
        error('observed.Y must have n*K rows.');
    end
    if isempty(Y)||size(Y,2)<1||any(~isfinite(Y(:)))
        error('observed.Y must contain at least one finite signal sample.');
    end
    d=sum(Y.^2,2)/size(Y,2);
    Sfull=max(d+d'-2*(Y*Y'/size(Y,2)),0);
    observed.S=Sfull;
    [P,observed]=prepare_problem(observed,o);
    return
end
if ~isequal(size(Sedge),[K,K,E]), error('S_edge must be K-by-K-by-E.'); end
if ~isequal(size(Rcopy),[K,K,n]), error('R_copy must be K-by-K-by-n.'); end
if any(~isfinite(Sedge(:)))||any(~isfinite(Rcopy(:)))
    error('S_edge and R_copy must be finite.');
end
for e=1:E
    Sedge(:,:,e)=max(0,0.5*(Sedge(:,:,e)+Sedge(:,:,e)'));
end
for i=1:n
    Ri=max(0,0.5*(Rcopy(:,:,i)+Rcopy(:,:,i)'));
    Ri(1:K+1:end)=0;
    Rcopy(:,:,i)=Ri;
end

if isfield(observed,'cost_scale') && ~isempty(observed.cost_scale)
    scale=observed.cost_scale;
elseif ~isempty(o.cost_scale)
    scale=o.cost_scale;
else
    vals=[Sedge(:);Rcopy(Rcopy>0)];
    vals=vals(isfinite(vals)&vals>0);
    scale=median(vals);
end
if ~isfinite(scale)||scale<=0, error('The cost scale must be positive.'); end
Sedge=Sedge/scale; Rcopy=Rcopy/scale;
P=struct('n',n,'K',K,'r',r,'E',E,'edge_i',ii,'edge_j',jj, ...
    'H',H,'S_edge',Sedge,'R_copy',Rcopy,'cost_scale',scale);
end

function [Theta,Gamma,B]=initialize(P,o,s)
if s==1 && ~isempty(o.initial_Theta), Theta=max(o.initial_Theta,0);
else, Theta=o.theta_scale*(0.25+rand(P.E,P.r)); end
if s==1 && ~isempty(o.initial_Gamma), Gamma=max(o.initial_Gamma,0);
else, Gamma=o.gamma_scale*(0.25+rand(P.n,P.r)); end
if s==1 && ~isempty(o.initial_B), B=project_simplex_columns(o.initial_B);
else
    B=random_interior_simplex(P.K,P.r,o.initial_uniform_mix);
end
end

function B=random_interior_simplex(K,r,uniform_mix)
% Independent Dirichlet(1) columns, mildly mixed with the barycenter.
% This is truth-independent and does not constrain later B iterates.
Z=-log(max(rand(K,r),realmin('double')));
B=Z./sum(Z,1);
B=(1-uniform_mix)*B+uniform_mix/K;
B=B./sum(B,1);
end

function G=reconstruct(Theta,Gamma,B,P)
Astate=zeros(P.n,P.n,P.r);
for m=1:P.r
    A=zeros(P.n);
    A(sub2ind([P.n,P.n],P.edge_i,P.edge_j))=Theta(:,m);
    A=A+A';
    Astate(:,:,m)=A;
end
Asupra=zeros(P.n*P.K);
for m=1:P.r
    bb=B(:,m)*B(:,m)';
    M=bb-diag(diag(bb));
    Asupra=Asupra+kron(bb,Astate(:,:,m))+kron(M,diag(Gamma(:,m)));
end
Asupra=0.5*(Asupra+Asupra');
Asupra(1:size(Asupra,1)+1:end)=0;
sa=sum(Theta,1)'; sc=sum(Gamma,1)';
view=B*diag(2*sa+sc)*B';
view=view-diag(diag(view));
G=struct('A_state',Astate,'A_supra',Asupra,'supra',Asupra, ...
    'view_graph',view,'state_cross_mass',sa,'state_copy_mass',sc);
end

function audit=gradient_check(P,o)
rng(o.seed+991,'twister');
Theta=.2+rand(P.E,P.r); Gamma=.2+rand(P.n,P.r);
B=.2+rand(P.K,P.r); B=B./sum(B,1);
[~,Gt,Gg,Gb]=objective_gradient(Theta,Gamma,B,P,o);
audit=struct();
[audit.theta_fd,audit.theta_analytic,audit.theta_relative_error]= ...
    directional_check(Theta,Gt,@(X)objective_only(X,Gamma,B,P,o),false);
[audit.gamma_fd,audit.gamma_analytic,audit.gamma_relative_error]= ...
    directional_check(Gamma,Gg,@(X)objective_only(Theta,X,B,P,o),false);
[audit.B_fd,audit.B_analytic,audit.B_relative_error]= ...
    directional_check(B,Gb,@(X)objective_only(Theta,Gamma,X,P,o),true);
audit.passed=max([audit.theta_relative_error,audit.gamma_relative_error, ...
    audit.B_relative_error])<o.gradient_check_tol;
end

function [fd,an,err]=directional_check(X,G,fun,tangent)
D=randn(size(X));
if tangent, D=D-mean(D,1); end
D=D/norm(D,'fro'); h=1e-6;
fd=(fun(X+h*D)-fun(X-h*D))/(2*h);
an=sum(G(:).*D(:));
err=abs(fd-an)/max([1,abs(fd),abs(an)]);
end

function r=projected_residual_nonnegative(X,G,L)
r=L*norm(X-max(X-G/L,0),'fro')/max(1,norm(X,'fro'));
end
function r=projected_residual_simplex(B,G,L)
r=L*norm(B-project_simplex_columns(B-G/L),'fro')/max(1,norm(B,'fro'));
end
function B=project_simplex_columns(Y)
B=zeros(size(Y));
for m=1:size(Y,2), B(:,m)=project_simplex(Y(:,m)); end
end
function x=project_simplex(y)
y=y(:); u=sort(y,'descend'); cssv=cumsum(u)-1;
rho=find(u-cssv./(1:numel(y))'>0,1,'last');
if isempty(rho), x=ones(size(y))/numel(y); return; end
x=max(y-cssv(rho)/rho,0);
x=x/sum(x);
end

function o=defaults(o)
d=struct('operation','solve','seed',1,'n_restarts',8,'max_iter',600, ...
 'eta_c',1,'alpha_a',0.01,'alpha_c',0.01,'beta_a',0.10,'beta_c',0.10, ...
 'tau_intra',0.10,'tau_inter',0.10,'lambda_B',0.01, ...
 'epsilon_intra',1e-5,'epsilon_inter',1e-5, ...
 'L_theta',1,'L_gamma',1,'L_B',1,'L_min',1e-10,'bt_factor',2, ...
 'sigma_theta',1e-6,'sigma_gamma',1e-6,'sigma_B',1e-6, ...
 'bt_relax',1.25,'max_backtracks',60,'bt_tol',0, ...
 'stationarity_tol',1e-6,'objective_tol',1e-9, ...
 'theta_scale',0.20,'gamma_scale',0.20,'initial_uniform_mix',0.10, ...
 'initial_Theta',[],'initial_Gamma',[],'initial_B',[], ...
 'cost_scale',[],'gradient_check_tol',5e-5,'verbose',false);
o=merge(d,o);
if o.beta_a<=0||o.beta_c<=0, error('beta_a and beta_c must be positive.'); end
epsnames={'epsilon_intra','epsilon_inter'};
for q=1:numel(epsnames)
    if o.(epsnames{q})<=0, error('%s must be positive.',epsnames{q}); end
end
nonneg={'eta_c','alpha_a','alpha_c','tau_intra','tau_inter','lambda_B'};
for q=1:numel(nonneg)
    if o.(nonneg{q})<0, error('%s must be nonnegative.',nonneg{q}); end
end
positive={'L_min','sigma_theta','sigma_gamma','sigma_B'};
for q=1:numel(positive)
    if o.(positive{q})<=0, error('%s must be positive.',positive{q}); end
end
if o.bt_tol~=0
    error('bt_tol must be zero for the stated sufficient-decrease guarantee.');
end
if o.initial_uniform_mix<0||o.initial_uniform_mix>=1
    error('initial_uniform_mix must lie in [0,1).');
end
end
function a=merge(a,b)
f=fieldnames(b);
for q=1:numel(f), a.(f{q})=b.(f{q}); end
end
