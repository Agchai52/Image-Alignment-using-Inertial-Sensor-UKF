function P = CovarianceHomo(h,covariance_x)
global PP1 PP2

PointsNum = size(PP1,2);
p1 = ones(3,PointsNum);
A  = zeros(2*PointsNum,9);

p1(1:2,:) = PP1;
H = reshape(h,[3,3])';
p2 = H*p1;
h = h/norm(h,2);

% p2 = PP2;
for i = 1:PointsNum
    A(2*i-1,:) = [p1(1,i) p1(2,i) 1 0 0 0 -p1(1,i)*p2(1,i) -p1(2,i)*p2(1,i) -p2(1,i)];
    A(2*i,:)   = [0 0 0 p1(1,i) p1(2,i) 1 -p1(1,i)*p2(2,i) -p1(2,i)*p2(2,i) -p2(2,i)];
end
AA = A'*A;
[V,D] = eig(AA);
J = 0;
for k = 2:9
    uk = V(:,k);
    lk = D(k,k);
    J = -(J+uk*uk'/lk);
end

S_small = 0;
S_biger = 0.5*(covariance_x(1,1)+covariance_x(2,2));
S = 0;
for i = 1:PointsNum
    a1 = A(2*i-1,:);
    a2 = A(2*i,:);
    foo = S_small*( h(1)^2 + h(2)^2 - 2*p2(1,i)*(h(1)*h(7)+h(2)*h(8)) )...
        + 2*S_biger*( p1(1,i)*h(7)*h(9) + p1(1,i)*p1(2,i)*h(7)*h(8) + p1(2,i)*h(8)*h(9) )...
        + ( S_small*p2(1,i)^2 + p1(1,i)^2*S_biger )*h(7)^2 ...
        + ( S_small*p2(1,i)^2 + p1(2,i)^2*S_biger )*h(8)^2 + S_biger*h(9)^2;
    fee = S_small*( h(4)^2 + h(5)^2 - 2*p2(2,i)*(h(4)*h(7)+h(5)*h(8)) )...
        + 2*S_biger*( p1(1,i)*h(7)*h(9) + p1(1,i)*p1(2,i)*h(7)*h(8) + p1(2,i)*h(8)*h(9) )...
        + ( S_small*p2(2,i)^2 + p1(1,i)^2*S_biger )*h(7)^2 ...
        + ( S_small*p2(2,i)^2 + p1(2,i)^2*S_biger )*h(8)^2 + S_biger*h(9)^2;
    feo = S_small*( (h(1)-p2(1,i)*h(7))*(h(4)-p2(2,i)*h(7)) + (h(2)-p2(1,i)*h(8))*(h(5)-p2(2,i)*h(8)) );
    foe = feo;
    
    S = S + a1'*a1*foo + a2'*a2*fee + a1'*a2*foe +a2'*a1*feo;
    
end

P = J*S*J;
P = 0.5*P+0.5*P';
P = P + 0.01*eye(9);
P = P/norm(P,2);
e = eig(P);
end