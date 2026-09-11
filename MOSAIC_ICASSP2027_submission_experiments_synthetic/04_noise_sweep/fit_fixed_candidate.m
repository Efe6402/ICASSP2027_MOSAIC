function est=fit_fixed_candidate(method,record,truth,row,cfg)
setup_noise_sweep_paths(cfg.root); method=upper(string(method));
switch method
    case "MOSAIC"
        names={'alpha_a','alpha_c','beta_a','beta_c','tau_intra','tau_inter', ...
            'lambda_B','eta_c','epsilon_intra','epsilon_inter','theta_scale','gamma_scale'};
        o=struct('verbose',false,'stationarity_tol',1e-5,'objective_tol',1e-10, ...
            'n_restarts',cfg.mosaic_restarts,'max_iter',cfg.mosaic_max_iter);
        for j=1:numel(names), o.(names{j})=double(row.(names{j})); end
        o.seed=double(mod(record.seed+1009*double(row.candidate_id),2^31-1));
        obs=struct('n',truth.n,'K',truth.K,'r',double(row.r),'S',record.S);
        t=tic; raw=solve_mosaic_crossview(obs,o); runtime=toc(t);
        est=common(method,raw.A_supra,runtime,raw.status,raw.objective);
        est.B=raw.B; est.Gamma=raw.Gamma; est.A_state=raw.A_state;
        est.iterations=raw.iterations;
        est.stationarity_max=raw.stationarity_max;
    case "PGL2021"
        p=struct('b1',double(row.beta1),'b2',double(row.beta2));
        t=tic; [Lp,Lq]=Learn_PGL(record.Sp,record.Sq,p); runtime=toc(t);
        Wp=symzero(max(0,diag(diag(Lp))-Lp)); Wq=symzero(max(0,diag(diag(Lq))-Lq));
        A=kron(eye(truth.K),Wp)+kron(Wq,eye(truth.n));
        est=common(method,A,runtime,'completed',NaN); est.Lp=Lp; est.Lq=Lq;
    case "ZW2024"
        o=struct('lambda',double(row.lambda),'rho_L',double(row.rho_L), ...
            'rho_C',double(row.rho_C),'max_iter',cfg.zw_max_iter, ...
            'tol',cfg.zw_tol,'backend',cfg.zw_backend);
        obs=struct('n',truth.n,'L',truth.K,'S',record.S);
        t=tic; raw=solve_zhang_wai_2024_paper(obs,o); runtime=toc(t);
        est=common(method,raw.A_effective,runtime,raw.status,raw.objective);
        est.A_primitive=raw.A_supra_primitive; est.iterations=raw.iterations;
    otherwise, error('Unknown method %s.',method);
end
est.candidate_id=double(row.candidate_id);
end
function e=common(method,A,runtime,status,obj)
e=struct('method',char(method),'A_primary',symzero(max(A,0)), ...
    'runtime',runtime,'status',char(status),'objective',obj);
end
function A=symzero(A), A=.5*(A+A'); A(1:size(A,1)+1:end)=0; end
