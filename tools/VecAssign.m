function [u , u_x , u_y] = VecAssign(v_x , v_y)

u_x = v_x;

u_y = v_y;

u = [u_x , u_y];
end