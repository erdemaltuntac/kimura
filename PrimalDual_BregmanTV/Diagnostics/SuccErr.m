function SuccErr_val = SuccErr(u_1 , u_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Erdem Altuntac
%
%                       Universite Libre de Bruxelles
%                       Department of Mathematics
%
%                       e-mail: erdem.altuntac@ulb.ac.be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_0 = u_0(:);
u_1 = u_1(:);

SuccErr_val = norm(u_1 - u_0);

end