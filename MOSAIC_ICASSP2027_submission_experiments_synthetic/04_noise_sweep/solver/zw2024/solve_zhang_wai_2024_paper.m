function result = solve_zhang_wai_2024_paper(observed,opts)
%SOLVE_ZHANG_WAI_2024_PAPER Paper-faithful AO solver for ICASSP 2024.
%
% Native learned objects:
%   A_1,...,A_L and C.
% Primitive supra-adjacency:
%   A_L + A_C, A_C=C kron I_n.
% High-order matched operator:
%   A_L + A_C + lambda(A_C A_L + A_L A_C).
%
% backend='cvx' reproduces the paper's stated numerical implementation.
% backend='exact_projection' solves the same convex graph subproblems in
% closed form and is included only as a fast mathematical diagnostic.

if nargin<2, opts=struct(); end
opts=defaults(opts);
n=observed.n; L=observed.L; D=n*L;
if n<2 || L<2, error('Require n>=2 and L>=2.'); end

if isfield(observed,'S') && isequal(size(observed.S),[D,D])
    S=0.5*(observed.S+observed.S');
elseif isfield(observed,'Y') && size(observed.Y,1)==D
    S=compute_zhang_wai_2024_distance_matrix(observed.Y);
else
    error('Provide raw Y or a D-by-D paper distance matrix S.');
end

if strcmpi(opts.backend,'cvx') && exist('cvx_begin','file')~=2
    error('CVX backend requested but CVX is unavailable on the MATLAB path.');
end

switch lower(char(opts.initialization))
    case 'uniform_feasible'
        C=uniform_graph(L);
        A=cell(L,1); for ell=1:L, A{ell}=uniform_graph(n); end
    otherwise
        error('Unknown initialization "%s".',opts.initialization);
end

obj_hist=nan(opts.max_iter+1,1);
change_hist=nan(opts.max_iter,3);
obj_hist(1)=objective(A,C,S,opts);
status='maximum_iterations';

for it=1:opts.max_iter
    ALold=blockdiag_cells(A); Cold=C;

    % Algorithm 1, Step 3.
    AC=kron(C,eye(n));
    P=S+opts.lambda*(AC*S+S*AC);
    for ell=1:L
        I=(ell-1)*n+(1:n);
        A{ell}=solve_graph_subproblem(P(I,I),opts.rho_L,n,opts.backend);
    end

    % Algorithm 1, Step 4.
    AL=blockdiag_cells(A);
    Q=S+opts.lambda*(AL*S+S*AL);
    G=zeros(L);
    for k=1:L
        Ik=(k-1)*n+(1:n);
        for ell=1:L
            Il=(ell-1)*n+(1:n);
            G(k,ell)=trace(Q(Ik,Il));
        end
    end
    C=solve_graph_subproblem(G,opts.rho_C,L,opts.backend);

    obj_hist(it+1)=objective(A,C,S,opts);
    change_AL=norm(AL-ALold,'fro');
    change_AC=sqrt(n)*norm(C-Cold,'fro'); % ||Delta(C kron I_n)||_F
    switch lower(char(opts.stopping_rule))
        case 'paper_absolute'
            ch=max(change_AL,change_AC);
        otherwise
            error('Paper reproduction supports stopping_rule="paper_absolute" only.');
    end
    change_hist(it,:)=[change_AL,change_AC,ch];
    if ch<=opts.tol
        status='converged';
        break;
    end
end

last=find(isfinite(obj_hist),1,'last');
AL=blockdiag_cells(A); AC=kron(C,eye(n));
primitive=AL+AC;
twohop=AC*AL+AL*AC;
effective=primitive+opts.lambda*twohop;

result=struct();
result.A_layers=A; result.C=C; result.A_L=AL; result.A_C=AC;
result.A_supra_primitive=0.5*(primitive+primitive');
result.A_twohop=0.5*(twohop+twohop');
result.A_effective=0.5*(effective+effective');
% IMPORTANT: Eq. (2) is the ordinary supra-adjacency, not Eq. (13).
result.A_supra=result.A_supra_primitive;
result.S=S; result.lambda=opts.lambda; result.rho_L=opts.rho_L; result.rho_C=opts.rho_C;
result.objective=obj_hist(last); result.objective_history=obj_hist(1:last);
result.change_history=change_hist(1:max(last-1,0),:);
result.iterations=last-1; result.status=status; result.backend=opts.backend;
result.options=opts;
end

function A=solve_graph_subproblem(D,rho,total_mass,backend)
D=0.5*(D+D'); d=size(D,1);
switch lower(char(backend))
    case 'cvx'
        cvx_begin quiet
            variable X(d,d) symmetric
            minimize( sum(sum(D.*X)) + rho*sum_square(X(:)) )
            subject to
                X >= 0;
                diag(X) == 0;
                sum(X(:)) == total_mass;
        cvx_end
        if ~(contains(lower(cvx_status),'solved'))
            error('CVX graph subproblem failed: %s',cvx_status);
        end
        A=full(X);
    case 'exact_projection'
        mask=triu(true(d),1);
        y=-D(mask)/(2*rho);
        x=project_scaled_simplex(y,total_mass/2);
        A=zeros(d); A(mask)=x; A=A+A';
    otherwise
        error('Unknown backend "%s".',backend);
end
A=0.5*(A+A'); A(1:d+1:end)=0;
end

function x=project_scaled_simplex(y,mass)
y=y(:);
if mass<0, error('Simplex mass must be nonnegative.'); end
if isempty(y), error('Graph has no strict-upper-triangle variables.'); end
u=sort(y,'descend'); cssv=cumsum(u)-mass;
rho=find(u-cssv./(1:numel(y))'>0,1,'last');
if isempty(rho), error('Scaled-simplex projection failed.'); end
theta=cssv(rho)/rho;
x=max(y-theta,0);
if abs(sum(x)-mass)>1e-9*max(1,mass)
    error('Scaled-simplex projection mass residual is too large.');
end
end

function value=objective(A,C,S,o)
n=size(A{1},1); AL=blockdiag_cells(A); AC=kron(C,eye(n));
Hhat=AL+AC+o.lambda*(AC*AL+AL*AC);
value=sum(sum(Hhat.*S));
for k=1:numel(A), value=value+o.rho_L*norm(A{k},'fro')^2; end
value=value+o.rho_C*norm(C,'fro')^2;
end

function A=uniform_graph(d)
A=(ones(d)-eye(d))/(d-1); % total mass exactly d
end

function AL=blockdiag_cells(A)
n=size(A{1},1); L=numel(A); AL=zeros(n*L);
for ell=1:L
    I=(ell-1)*n+(1:n); AL(I,I)=A{ell};
end
end

function o=defaults(o)
d=struct('lambda',0.1,'rho_L',1,'rho_C',15,'max_iter',100, ...
    'tol',1e-3,'stopping_rule','paper_absolute','backend','cvx', ...
    'initialization','uniform_feasible','store_objective_history',true);
f=fieldnames(d); for k=1:numel(f), if ~isfield(o,f{k}), o.(f{k})=d.(f{k}); end, end
if o.lambda<0 || o.rho_L<=0 || o.rho_C<=0, error('Invalid lambda/rho values.'); end
end
