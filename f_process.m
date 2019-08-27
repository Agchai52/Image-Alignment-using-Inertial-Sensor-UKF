function Xaug_next = f_process(Xaug)
    global dt
    h = dt;
    Xw = Xaug(1:6,1);
    F1 = h*fw(Xw);
    F2 = h*fw(Xw+F1/2);
    F3 = h*fw(Xw+F2/2);
    F4 = h*fw(Xw+F3);
    Xw_next = Xw + 1/6*(F1+2*F2+2*F3+F4);
    clear F1 F2 F3 F4
    
    R  = convertE2R(Xw_next(4:6,1));
    Xa = Xaug(7:15,1);
    Xa(1:3,1) = inv(R)*Xa(1:3);
    F1 = h*fa(Xa);
    F2 = h*fa(Xa+F1/2);
    F3 = h*fa(Xa+F2/2);
    F4 = h*fa(Xa+F3);
    Xa_next = Xa + 1/6*(F1+2*F2+2*F3+F4);
    
    Xaug_next = [Xw_next;Xa_next];
    
end

function X_dot = fw(X0)  
   %X = [wx wy wz anglex angley anglez]
    A = [0 0 0 0 0 0;
         0 0 0 0 0 0;
         0 0 0 0 0 0;
         1 0 0 0 0 0; % dtheta_x/dt = omega_x
         0 1 0 0 0 0; % dtheta_y/dt = omega_y
         0 0 1 0 0 0]; % dtheta_z/dt = omega_z
     X_dot = A*X0;        
end

function X_dot = fa(X0)  
   %X = [ax ay az vx vy vz tx ty tz];
    A = [0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 0 0 0 0;
         1 0 0 0 0 0 0 0 0; % dv_x/dt = a_x
         0 1 0 0 0 0 0 0 0; 
         0 0 1 0 0 0 0 0 0;
         0 0 0 1 0 0 0 0 0; % dt_x/dt = v_x
         0 0 0 0 1 0 0 0 0;
         0 0 0 0 0 1 0 0 0]; 
     X_dot = A*X0;        
end

function R = convertE2R(eulerangle)
   Rx = [1 0 0;
         0 cos(eulerangle(1)) -sin(eulerangle(1));
         0 sin(eulerangle(1))  cos(eulerangle(1))];
   Ry = [ cos(eulerangle(2)) 0 sin(eulerangle(2));
          0 1 0;
         -sin(eulerangle(2)) 0   cos(eulerangle(2))];
   Rz = [ cos(eulerangle(3)) -sin(eulerangle(3)) 0;
          sin(eulerangle(3))  cos(eulerangle(3)) 0
          0 0 1];
   R = Rz*Ry*Rx;              
end