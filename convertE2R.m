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