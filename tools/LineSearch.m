function [u_rec] = LineSearch(u0, u1, lambda)
%%                      LineSearch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Erdem Altuntac
%
%                       Universite Libre de Bruxelles
%                       Department of Mathematics
%
%                       e-mail: erdem.altuntac@ulb.ac.be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u0 = u0(:);
u1 = u1(:);

u_rec = (1 - lambda).*u0 + lambda .* u1;

end