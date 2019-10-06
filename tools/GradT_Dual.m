function [DT_w, DTw_x , DTw_y] = GradT_Dual(w_x , w_y , DxT , DyT)
%%                       GradT_Dual
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
% Subdifferential of non-smooth total variation is formulated by the 
% transpose of the gradient matrix applied on the some dual variable.
%
% INPUT:
% -------------------------------------------------------------------------
% w                  :       Dual variable
% DxT , DyT , DzT    :       Tranpose of the derivatives in all directions
% 
% OUTPUT:
% -------------------------------------------------------------------------
% DT_w    :   Transpose of the gradient applied on the dual variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [m,n] = size(DxT);

DTw_x = DxT * w_x;
DTw_y = DyT * w_y;


DT_w = DTw_x + DTw_y;%[DTw_x , DTw_y , DTw_z];

% if(nargin == 3)
%     DzT = zeros(m,n);
%     DTw_x = DxT * w_x;
%     DTw_y = DyT * w_y;
%     DTw_z = DzT * w_z;
%     
%     DT_w = [DTw_x , DTw_y , DTw_z];
%     
%     %DT_w = DxT * w_x + DyT * w_y + DzT * w_z;
% else
%         
%     DT_w = [DTw_x , DTw_y , DTw_z];
% end

end