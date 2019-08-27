function [R1,T1,t1,angle1,R2,T2,t2,angle2] = decompHomoMatrix2(H)

% H = [0.98196095,0.16305955,-6.1110927e-10;-0.17240666,0.99494064,-2.2870511e-11;45.054413,0.53569812,0.99999994]

% theta = 0.5;
% Rx = eye(3);%[1 0 0; 0 cosd(theta) -sind(theta);0  sind(theta) cosd(theta)]
% Ry = [cosd(theta) 0 sind(theta) ;0 1 0;-sind(theta) 0 cosd(theta)];
% Rz = [cosd(theta) -sind(theta) 0;sind(theta) cosd(theta) 0; 0 0 1];
% Rall = Rx*Ry*Rz
% H = Rall*[1 0 0;0 1 0;5 5 1]
%H = [0.998782025129912,-0.0348782368720627,0.0348994967025010;0.0360954697847683,0.998739518419950,-0.0348782368720627;2.07221144449168,2.96255755130049,0.963946307918726]
Rp = cell(4,1);
tp = cell(4,1);
R = cell(4,1);
T = cell(4,1);

[U,S,V] = svd(H);


d1 = S(1,1);
d2 = S(2,2);
d3 = S(3,3);
s = det(U)*det(V);

d = s*d2;

x1 = sqrt((d1^2-d2^2)/(d1^2-d3^2));
x3 = sqrt((d2^2-d3^2)/(d1^2-d3^2));

np1 = [x1;0;x3];
np2 = [x1;0;-x3];
np3 = [-x1;0;x3];
np4 = [-x1;0;-x3];

np = [np1,np2,np3,np4];

n = [V*np1,V*np2,V*np3,V*np4];
for i=1:4
    if sum(n(:,i))*d>0
        sin_theta = sign(np(1,i)*np(3,i))*sqrt((d1^2-d2^2)*(d2^2-d3^2))/((d1+d3)*d2);
        cos_theta = (d2^2+d1*d3)/((d1+d3)*d2);
        Rp{i} = [cos_theta 0 -sin_theta;0 1 0;sin_theta 0 cos_theta];
        tp{i} = (d1-d3)*np(:,i);
        %2*i-1

    elseif sum(n(:,i))*d<0
        sin_theta = sign(np(1,i)*np(3,i))*sqrt((d1^2-d2^2)*(d2^2-d3^2))/((d1-d3)*d2);
        cos_theta = (d1*d3-d2^2)/((d1-d3)*d2);
        Rp{i} = [cos_theta 0 sin_theta;0 -1 0;sin_theta 0 -cos_theta];
        tp{i} = (d1+d3)*np(:,i);
        %2*i
    end

end
j=0; 
for i = 1:4
    R{i} = s*U*Rp{i}*V';
    T{i} = U*tp{i};
    temp = -inv(R{i})*T{i}*n(:,i)';
    if n(3,i)>0 && j==0 %
%         fprintf('The %d th result is a possible solution\n',i);
        R1 = R{i};
        T1 = T{i};
        t1 = temp(3,:)';
        n1 = n(:,i);
        angle1 = convertR2E(R1);
        j = j+1;
    elseif n(3,i)>0 && j==1 
%         fprintf('The %d th result is a possible solution\n',i);
        R2 = R{i};
        T2 = T{i};
        t2 = temp(3,:)';
        n2 = n(:,i);
        angle2 = convertR2E(R2);
    end
end
% R1 = R1/R1(3,3);
% R2 = R2/R2(3,3);

% H1 = R1-T1*n1';
% H1 = R1*(eye(3)+[0;0;1]*t1');
% lambda1 = H1(3,3);
% %   H1 = H1/H1(3,3)
% %H2 = R2-T2*n2';
% H2 = R2*(eye(3)+[0;0;1]*t2');
% %   H2 = H2/H2(3,3)
% lambda2 = H2(3,3); 


% eigH = [d1 d2 d3]';
end

function eulerangle = convertR2E(R)
   eulerangle = zeros(3,1);
   eulerangle(1,1) = atan(R(3,2)/R(3,3)); % euler angle for x-axis
   eulerangle(2,1) = asin(-R(3,1)); % euler angle for y-axis
   eulerangle(3,1) = atan(R(2,1)/R(1,1)); % euler angle for z-axis
end