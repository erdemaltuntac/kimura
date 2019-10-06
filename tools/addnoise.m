function [v_delta , delta] = addnoise(v_true)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                        Erdem Altuntac
%                        Universite Libre de Bruxelles
%                        Department of Mathematics
%
%                        Erdem.Altuntac@ulb.ac.be
%                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%v_std = std(v_true(:));
%delta = .05;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Additive Noise %%%%%%%%%%
%v_delta = v_true + delta .* rand(length(v_true),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Multiplicative %%%%%%%%%%
%v_delta = delta .* rand(length(v_true),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_v = length(v_true);
err_lev   = 200;
delta     = err_lev/100 * norm(v_true(:)) / sqrt(L_v);
v_delta   = v_true(:) + delta .* rand(length(v_true),1);

end