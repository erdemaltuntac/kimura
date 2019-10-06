function [Mfit , Mfit_norm] = Misfit(A,x,b)

x = x(:);
b = b(:);

Mfit = A * x - b;
Mfit_norm = norm(Mfit);

end