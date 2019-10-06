function [w_new , w_new_x , w_new_y] = DualUpdate(w_x , w_y , nu , regpar, u_x , u_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Erdem Altuntac
%
%                       Universite Libre de Bruxelles
%                       Department of Mathematics
%
%                       e-mail: erdem.altuntac@ulb.ac.be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dual update is meant to be global for multiple algorithms.
% In case of BregmanTV algorithm, (u_x,u_y,u_z) is the gradient of u in
% each direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_x = w_x + nu .* u_x;
v_y = w_y + nu .* u_y;

[DualVec , DualVec_x , DualVec_y] = VecAssign(v_x , v_y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Prox Input Preparation %%%%%%%%%%%%%

PrxInpt.nu = nu;
PrxInpt.regpar = regpar;
PrxInpt.DualVec = DualVec;
PrxInpt.DualVec_x = DualVec_x;
PrxInpt.DualVec_y = DualVec_y;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[w_new , w_new_x , w_new_y] = prox_l1_conj(PrxInpt);


end