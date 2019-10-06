function [Grad , Grad_T, Dx , Dy] = FiniteDiff2D(Geometry)
%%                       FiniteDiff2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Erdem Altuntac
%
%                       Universite Libre de Bruxelles
%                       Department of Mathematics
%
%                       e-mail: erdem.altuntac@ulb.ac.be
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Two dimensional finite difference method for gradient with respect to 
% the forward difference approach.
%
%   INPUT:
%   Geometry : Area of interest - domain
% -------------------------------------------------------------------------
%   OUTPUT:
%   Dx , Dy  : Difference matrices for each spatial coordinate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[Params] = ProblemParams();
% Lx = Geometry.Lx;
% Ly = Geometry.Ly;
nx = Geometry.Nx; 
ny = Geometry.Ny;


% x0 = Geometry.x0;
% y0 = Geometry.y0;

hx = Geometry.DeltaX;%(Lx - x0) / (nx - 1); % Mesh length
hy = Geometry.DeltaY;%(Ly - y0) / (ny - 1);
%% Lexicographical scheme w.r.t. central finite diff

% x_coeffs = [-ones(nx,1) , ones(nx,1)];
% y_coeffs = [-ones(ny,1) , ones(ny,1)];
% 
% x_diff = (1/(hx)).*spdiags(x_coeffs, [0 1], nx - 1 , nx);
% y_diff = (1/(hy)).*spdiags(y_coeffs, [0 1], ny - 1 , ny);

x_coeffs = [-ones(nx,1) , zeros(nx,1) , ones(nx,1)];
y_coeffs = [-ones(ny,1) , zeros(ny,1) , ones(ny,1)];
% z_coeffs = [-ones(nz,1) , zeros(nz,1) , ones(nz,1)];
% 
x_diff = (1/(2*hx)).*spdiags(x_coeffs, [0 1 2], nx - 2 , nx);
y_diff = (1/(2*hy)).*spdiags(y_coeffs, [0 1 2], ny - 2 , ny);
% z_diff = (1/(2*hz)).*spdiags(z_coeffs, [0 1 2], nz - 2 , nz);

% Difference matrices

Dx = kron(speye(nx) , x_diff);
Dy = kron(y_diff , speye(ny));

Grad = [Dx , Dy];

Grad_T = -Grad';


end