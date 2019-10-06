function [Grad_u , u_x , u_y , TV_u] = DiscrDiff_Val(u, Dx , Dy)

%% Gradient of the given vector in terms of finite difference

u_x = Dx * u(:);
u_y = Dy * u(:);

Grad_u = [u_x , u_y];%u_x + u_y + u_z;

%% TV of the function is the l^1 of Euclidean norm of the gradient
Grad_u_Euclid = sqrt(u_x.^2 + u_y.^2);
TV_u = norm(Grad_u_Euclid , 1);

end