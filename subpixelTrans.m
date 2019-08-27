function [u,v] = subpixelTrans(Referece_L,Alternative_L,u_hat,v_hat,n)
% Subpixel accurate Translation
%--input
% Referece_L; reference for T
% Alternative_L; alternative for I
% u_hat = 0; % Integral alignment detected 
% v_hat = 0; % Integral alignment detected 
% n  = 8; % size of tile according to supplement
%--output
% [u,v]' = [u_hat,v_hat]'+ mu   % Equ.(13) - Equ. (30) in supplement
% u,v  % Subpixel alignment



if abs(u_hat)>=1 && abs(u_hat)>=1

[Lm,Ln]=size(Referece_L);
dm = fix(Lm/4);
dn = fix(Ln/4);
Lm = fix(Lm/2);
Ln = fix(Ln/2);
Ty = Lm-dm;%[Lm-dm,Lm,Lm+dm]; % location for tiles
Tx = Ln-dn;%[Ln-dn,Ln,Ln+dn];

du = [-1 0 1]; % possible offset
dv = [-1 0 1];

u_all = zeros(length(Tx),length(Ty));
v_all = zeros(length(Tx),length(Ty));

I_uv = zeros(length(du),length(dv)) ;

Ref = Referece_L;
Alt = Alternative_L;

D2_sub = zeros(length(du),length(dv)); % distance

for i = 1:length(Ty)
    for j = 1:length(Tx)
        dx = Tx(i);
        dy = Ty(j);
        
        T = Ref((1+dy):(n+dy),(1+dx):(n+dx));
        D2term1 = sum(sum(T.*T)) * ones(length(du),length(du));
        
        I_all = Alt((1+dy-1+u_hat):(n+dy+1+u_hat),(1+dx-1+v_hat):(n+dx+1+v_hat));
        h = ones(n,n);
        J = imfilter(I_all.*I_all,h);
        D2term2 = J(n/2:n/2-1+length(du),n/2:n/2-1+length(dv));
        
        m_iall = n+2;
        C = real(ifft2(fft2(I_all).* fft2(rot90(T,2),m_iall,m_iall)));
        D2term3 = C(n:n-1+length(du),n:n-1+length(dv));
        
        D2_sub = D2term1 + D2term2 - 2*D2term3;
        

%% Figure
%         dy = 270;
%         dx = 515;
%         
%         T2 = Ref((1+dy):(n+dy),(1+dx):(n+dx));
%         D2term1f = sum(sum(T2.*T2)) * ones(2*n+1,2*n+1);
%         
%         I_allf = Alt((1+dy-n+u_hat):(n+dy+n+u_hat),(1+dx-n+v_hat):(n+dx+n+v_hat));
%         h = ones(n,n);
%         J = imfilter(I_allf.*I_allf,h);
%         D2term2f = J(n/2:n/2-1+2*n+1,n/2:n/2-1+2*n+1);
%         
%         m_iall = n+2*n;
%         C = real(ifft2(fft2(I_allf).* fft2(rot90(T2,2),m_iall,m_iall)));
%         D2term3f = C(n:n-1+2*n+1,n:n-1+2*n+1);
%         
%         D2 = D2term1f + D2term2f - 2*D2term3f;
%         
%         figure; imshow(T2/255);
%         figure; imshow(I_allf/255);
%         figure; contourf(1:(2*n+1),1:(2*n+1),D2)

%% subpixelTrans


        FA11 = [1 -2 1;2 -4 2;1 -2 1]/4;
        FA22 = [1 2 1;-2 -4 -2;1 2 1]/4;
        FA12 = [1 0 -1;0 0 0;-1 0 1]/4;
        FA21 = FA12;

        Fb1 = [-1 0 1;-2 0 2;-1 0 1]/8;
        Fb2 = [-1 -2 -1;0 0 0;1 2 1]/8;

        Fc = [-1 2 -1;2 12 2;-1 2 -1]/16;

        A = [sum(dot(FA11,D2_sub)) sum(dot(FA12,D2_sub));sum(dot(FA21,D2_sub)) sum(dot(FA22,D2_sub))];
        b = [sum(dot(Fb1,D2_sub)); sum(dot(Fb2,D2_sub))];
        c = dot(Fc,D2_sub);

        A(1,1) = max(0,A(1,1));
        A(2,2) = max(0,A(2,2));


        det_A = A(1,1)*A(2,2)-A(1,2)^2;
        if det_A<0
            A(1,2)=0;
            A(2,1)=0;
        end

        mu = -inv(A)*b;
        
        if abs(mu(1))>=1 
            u_all(i) = u_hat;            
        else
            u_all(i) = mu(1)+u_hat;            
        end
        
        if abs(mu(2))>=1
            v_all(j) = v_hat;
        else
            v_all(i) = mu(2)+v_hat;
        end
    end
end

u = sum(sum(u_all))/length(Tx)/length(Ty);
v = sum(sum(v_all))/length(Tx)/length(Ty);

else
    u = u_hat;
    v = v_hat;
end