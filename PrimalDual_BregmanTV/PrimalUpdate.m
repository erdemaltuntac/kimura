function u_new = PrimalUpdate(u_old , Grad_Mfit , regpar , mu , w_dual)

% We pass the input onto the proximal mapping of the indicator function

Backward = Grad_Mfit + regpar .* w_dual;
Forward = u_old(:) - mu .* Backward;

u_new = prox_indicator(Forward);

  

end


