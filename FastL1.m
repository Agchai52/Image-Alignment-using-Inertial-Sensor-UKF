function [u,v] = FastL1(Referece_L,Alternative1_L,u0,v0,n)
% u0 = u1;
% v0 = v1;
% n  = 16;
% u_hat = 1;
% v_hat = 1;

[Lm,Ln]=size(Referece_L);
Tx = [3*n:n:Lm-3*n];
Ty = [3*n:n:Ln-3*n];

u_all = zeros(length(Tx),length(Ty));
v_all = zeros(length(Tx),length(Ty));


T = Referece_L;
I = Alternative1_L;

du = -1:1;
dv = -1:1;
D2 = zeros(length(du),length(dv));

for i = 1:length(Tx)
    for j = 1:length(Ty)
        dx = Tx(i);
        dy = Ty(j);
    try
    for u_hat = du
        for v_hat = dv
            for x=(1+dx):(n+dx)
                for y=(1+dy):(n+dy)
                    D2(u_hat+du(end)+1,v_hat+du(end)+1) = D2(u_hat+du(end)+1,v_hat+du(end)+1)+abs(T(x,y) - I(x+u_hat+u0,y+v_hat+v0));    
                end
            end
        end
    end
    [row,column]=find(D2==min(min(D2)));
    
    u_hat = du(row);
    v_hat = dv(column);

    u_all(i,j)  = u0+u_hat;
    v_all(i,j)  = v0+v_hat;

%     if abs(u_all(i,j))>=5 && abs(v_all(i,j))<5
%         u_all(i,j) = u0;
%     elseif abs(u_all(i,j))<5 && abs(v_all(i,j))>=5
%         v_all(i,j) = v0;
%     elseif abs(u_all(i,j))>=5 && abs(v_all(i,j))>=5
%         u_all(i,j) = u0;
%         v_all(i,j) = v0;
%     end
        
    catch
        u_all(i,j) = u0;
        v_all(i,j) = v0;
    continue;
    end
    end
end
u = round(sum(sum(u_all))/length(Tx)/length(Ty));
v = round(sum(sum(v_all))/length(Tx)/length(Ty));




end


