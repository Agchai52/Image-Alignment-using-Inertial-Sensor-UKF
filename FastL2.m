function [u,v] = FastL2(Referece_L,Alternative_L,u0,v0,n,level)
% Fast L2 residual computation
%--input
% Referece_L; reference for T
% Alternative_L; alternative for I
% u0 = 0; % Initial alignment inherited by the tile from coarser level
% v0 = 0; %
% n  = 8; % size of tile according to supplement
%--output
% u min of D2 = T_norm + I_norm - 2*ifft2(fft2(I).* fft2(T)); Equ.12 sup
% v

[Lm,Ln]=size(Referece_L);

if level==3
    Ty = fix(Lm/2)-4:2:fix(Lm/2)+4;
    Tx = fix(Ln/2)-4:2:fix(Ln/2)+4;
else
    Ty = [3*n:n:Lm-4*n]; % location for tiles
    Tx = [3*n:n:Ln-4*n];
end

% if Ty==[] || Tx==[]
%     Ty = fix(Lm/2)-4;
%     Tx = fix(Ln/2)-4;
% end

du = -4:4; % possible offset
dv = -4:4;


u_all = zeros(length(Tx),length(Ty));
v_all = zeros(length(Tx),length(Ty));

Ref = Referece_L;
Alt = Alternative_L;

D2 = zeros(length(du),length(dv)); % distance

for i = 1:length(Ty)
    for j = 1:length(Tx)
        dx = Tx(j);
        dy = Ty(i);
        
        
        T = Ref((1+dy):(n+dy),(1+dx):(n+dx));
        D2term1 = sum(sum(T.*T)) * ones(length(du),length(du));
        
        I_all = Alt((1+dy-4+u0):(n+dy+4+u0),(1+dx-4+v0):(n+dx+4+v0));
        h = ones(n,n);
        J = imfilter(I_all.*I_all,h);
        D2term2 = J(n/2:n/2-1+length(du),n/2:n/2-1+length(dv));
        
        m_iall = n+8;
        C = real(ifft2(fft2(I_all).* fft2(rot90(T,2),m_iall,m_iall)));
        D2term3 = C(n:n-1+length(du),n:n-1+length(dv));
        
        D2 = D2term1 + D2term2 - 2*D2term3;
        
        [row,column]=find(D2==min(min(D2)));
        
        if max(size(row))>1 || max(size(column))>1
            u_all(i,j) = du(min(row))+u0;
            v_all(i,j) = dv(min(column))+u0;
        else
            u_all(i,j) = du(row)+u0;
            v_all(i,j) = dv(column)+v0;
        end
       
    end
end

u = round(sum(sum(u_all))/length(Tx)/length(Ty));
v = round(sum(sum(v_all))/length(Tx)/length(Ty));

end