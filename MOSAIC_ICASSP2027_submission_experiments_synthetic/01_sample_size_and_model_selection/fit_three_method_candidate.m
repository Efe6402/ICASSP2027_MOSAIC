function est=fit_three_method_candidate(method,record,truth,row,cfg,stage)
%FIT_THREE_METHOD_CANDIDATE Fit one candidate to one signal realization.
setup_three_method_paths(); method=upper(string(method));
switch method
    case "MOSAIC"
        names={'alpha_a','alpha_c','beta_a','beta_c','tau_intra', ...
            'tau_inter','lambda_B','eta_c','epsilon_intra', ...
            'epsilon_inter','theta_scale','gamma_scale'};
        o=struct('verbose',false,'stationarity_tol',1e-5,'objective_tol',1e-10);
        for j=1:numel(names), o.(names{j})=double(row.(names{j})); end
        if strcmpi(stage,'screen')
            o.n_restarts=cfg.mosaic_screen_restarts; o.max_iter=cfg.mosaic_screen_max_iter;
        else
            o.n_restarts=cfg.mosaic_refine_restarts; o.max_iter=cfg.mosaic_refine_max_iter;
        end
        o.seed=double(mod(record.seed+1009*double(row.candidate_id),2^31-1));
        obs=struct('n',truth.n,'K',truth.K,'r',double(row.r),'S',record.S);
        t=tic; raw=solve_mosaic_crossview(obs,o); runtime=toc(t);
        est=common(method,raw.A_supra,runtime,raw.status,raw.objective);
        est.B=raw.B; est.Gamma=raw.Gamma; est.A_state=raw.A_state;
        est.raw=raw; est.A_primitive=raw.A_supra;
    case "PGL2021"
        param=struct('b1',double(row.beta1),'b2',double(row.beta2));
        t=tic; [Lp,Lq]=Learn_PGL(record.Sp,record.Sq,param); runtime=toc(t);
        Wp=max(0,diag(diag(Lp))-Lp); Wp=symzero(Wp);
        Wq=max(0,diag(diag(Lq))-Lq); Wq=symzero(Wq);
        A=kron(eye(truth.K),Wp)+kron(Wq,eye(truth.n));
        est=common(method,A,runtime,'completed',NaN);
        est.W_physical=Wp; est.W_view=Wq; est.Lp=Lp; est.Lq=Lq;
        est.A_primitive=est.A_primary;
    case "ZW2024"
        o=struct('lambda',double(row.lambda),'rho_L',double(row.rho_L), ...
            'rho_C',double(row.rho_C),'max_iter',cfg.zw_max_iter, ...
            'tol',cfg.zw_tol,'backend',cfg.zw_backend);
        obs=struct('n',truth.n,'L',truth.K,'S',record.S);
        t=tic; raw=solve_zhang_wai_2024_paper(obs,o); runtime=toc(t);
        est=common(method,raw.A_effective,runtime,raw.status,raw.objective);
        est.A_primitive=raw.A_supra_primitive; est.A_effective=raw.A_effective;
        est.A_layers=raw.A_layers; est.C=raw.C; est.raw=raw;
    otherwise
        error('Unknown method %s.',method);
end
est.candidate_id=double(row.candidate_id);
end
function est=common(method,A,runtime,status,obj)
est=struct('method',char(method),'A_primary',symzero(max(A,0)), ...
    'runtime',runtime,'status',char(status),'objective',obj);
end
function A=symzero(A)
A=0.5*(A+A'); A(1:size(A,1)+1:end)=0;
end
