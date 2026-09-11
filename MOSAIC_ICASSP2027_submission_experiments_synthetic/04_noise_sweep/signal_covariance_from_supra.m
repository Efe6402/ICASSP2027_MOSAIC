function [SigmaClean,Sigma]=signal_covariance_from_supra(A,c)
%SIGNAL_COVARIANCE_FROM_SUPRA Match the default heat/filter generator.
d=sum(A,2); L=diag(d)-A; lm=max(real(eig(L))); if lm>0, L=L/lm; end
[V,D]=eig(.5*(L+L')); ell=max(real(diag(D)),0);
if strcmpi(c.signal_model,'heat_diffusion')
    h=exp(-c.heat_time*ell);
else
    h=(1+c.filter_strength*ell).^(-c.filter_order);
end
H=V*diag(h)*V'; SigmaClean=c.signal_variance*(H*H');
SigmaClean=.5*(SigmaClean+SigmaClean');
Sigma=SigmaClean+c.observation_noise_variance*eye(size(A,1));
Sigma=.5*(Sigma+Sigma');
end
