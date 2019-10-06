function alpha_1 = RegPar(Mfit_norm, tau, u, delta, k)

% q_sqr = (q_rough)^2;
% alpha_1 = alpha_0 * q_sqr^(k-1);
% a = exp(-DiscRad);
% alpha_1 = a;%min(a , DiscRad^2);

yy = IndexFunc(Mfit_norm);

u_mag = norm(u);

tt = tau^2 - 1;

st = 1*((tt)*yy - .5*u_mag);

st_inv = 1/st;

alpha_0 = st_inv*Mfit_norm^2;

alpha_1 = delta^2 * alpha_0;

yy_dbl = 2*yy;
 
alpha_1 = max((delta^2/(yy_dbl))*(tt),min(alpha_1 , delta^2/yy));
 
end



