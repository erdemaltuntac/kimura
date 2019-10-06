function [w_new , w_new_x , w_new_y] = prox_l1_conj(PrxInpt)
%%                      prox_l1_conj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Erdem Altuntac
%
%                       Universite Libre de Bruxelles
%                       Department of Mathematics
%
%                       e-mail: erdem.altuntac@ulb.ac.be
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prox operator of the l^1-norm conjugation; projection over the interval
% [-nu , nu] where nu is given with the struct PrxInput. The projection is
% defined as follows, after the extraction DualVec = PrxInpt.DualVec,
%                  -
%                 / nu,            for (DualVec)_i > nu
%   (w_new)_i =  <  DualVec_i,     for abs(DualVec)_i <= nu 
%                 \ -nu,           for (DualVec)_i < nu
%                  -
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract the necessary input

w_x = PrxInpt.DualVec_x;
w_y = PrxInpt.DualVec_y;

regpar = PrxInpt.regpar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_norm = sqrt(w_x.^2 + w_y.^2);
[w0_new , w0_new_x , w0_new_y] = VecAssign(w_x , w_y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (w_norm > 1)
    
    w1_new_x = 1/w_norm * w0_new_x;
    w1_new_y = 1/w_norm * w0_new_y;
    
    [w_new , w_new_x , w_new_y] = VecAssign(w1_new_x , w1_new_y);
else
    
    [w_new , w_new_x , w_new_y] = VecAssign(w_x , w_y);
    
end


end
