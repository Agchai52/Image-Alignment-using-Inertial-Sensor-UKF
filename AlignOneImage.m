function [ img_aligned,du,dv ] = AlignOneImage(img_reference, img_offset)

[m,n]=size(img_reference);
%% Gaussian Pyramid
w=fspecial('gaussian',[2 2]);

Referece_L0=img_reference; %L0 means level 0 the bottom of Pyramid
Referece_L1=imresize(imfilter(Referece_L0,w),[m/2 n/2]);
Referece_L2=imresize(imfilter(Referece_L1,w),[m/2/4 n/2/4]);
Referece_L3=imresize(imfilter(Referece_L2,w),[m/2/4/4 n/2/4/4]);

Alternative1_L0=img_offset;
Alternative1_L1=imresize(imfilter(Alternative1_L0,w),[m/2 n/2]);
Alternative1_L2=imresize(imfilter(Alternative1_L1,w),[m/2/4 n/2/4]);
Alternative1_L3=imresize(imfilter(Alternative1_L2,w),[m/2/4/4 n/2/4/4]);

%% Fast subpixel L2 alignment
% Alternative1
u0 = 0;
v0 = 0;

%--Level 3: Fast subpixel L2 alignment
[u3,v3] = FastL2(Referece_L3,Alternative1_L3,u0,v0,8,3);
%--Level 2: Fast subpixel L2 alignment
[u2,v2] = FastL2(Referece_L2,Alternative1_L2,u3*4,v3*4,16,2);
%--Level 1: Fast subpixel L2 alignment
[u1,v1] = FastL2(Referece_L1,Alternative1_L1,u2*4,v2*4,16,1);
%--Level 0: L1 alignment
[u_hat,v_hat] = FastL1(Referece_L0,Alternative1_L0,u1*2,v1*2,16);

% % % % Subpixel accurate translation 
% % % % [u,v] = subpixelTrans(Referece_L0,Alternative1_L0,u_hat,v_hat,16);
% % % %fprintf('Estimated offset is (dy,dx)=(%d,%d)\n',u,v);
% % % 
% % % %% Cancel offset
% % % % I1E = MatrixExpand(img_offset);
% % % % [me,ne] = size(I1E);
% % % % 
% % % % du = round(u*4);
% % % % dv = round(v*4);
% % % 
% % % % if du>=0 && dv>=0 
% % % %     I1E = I1E(du+1:me,dv+1:ne);
% % % % elseif du<0 && dv>=0
% % % %     I1E = I1E(1:me+du,dv+1:ne);
% % % % elseif du>=0 && dv<0
% % % %     I1E = I1E(du+1:me,1:ne+dv);
% % % % elseif du<0 && dv<0
% % % %     I1E = I1E(1:me+du,1:ne+dv);
% % % % end
% % % 
% % % %% Downsampling 
% % % % tmp1=downsample(I1E,4);
% % % % tmp2=permute(tmp1,[2,1]);
% % % % tmp3=downsample(tmp2,4);
% % % % img_aligned = permute(tmp3,[2,1]);

%% No subpixel
du = u_hat;
dv = v_hat;


if du>=0 && dv>=0 
    I1E = img_offset(du+1:m,dv+1:n);
elseif du<0 && dv>=0
    I1E = img_offset(1:m+du,dv+1:n);
elseif du>=0 && dv<0
    I1E = img_offset(du+1:m,1:n+dv);
elseif du<0 && dv<0
    I1E = img_offset(1:m+du,1:n+dv);
end

if du~=0 || dv~=0
    I1E(m,n)=0;
end
img_aligned = I1E;
end

