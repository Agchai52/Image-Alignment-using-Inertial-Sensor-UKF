clc
clear all
close all

%% input
H0 = [1.000721311073639   0.009702829197208   0.000001094268079;
     -0.012149941577035   0.996074783341494  -0.000002314204994;
     12.613398039458495 -11.191775954098418   1.000000000000000]

 T2 =[1.0057312   0.0327093   0.0000079;
     -0.0458396   1.0013561  -0.0000114;
     24.0083961 -29.1256638   1.0000000] 
% H0 = T2;
 Point2 =[716.11823,841.70392;372.89124,636.93542;543.99042,746.05664;521.09406,753.20056;501.56699,740.99310;427.90366,640.53552;528.47461,707.41132;405.96606,825.64563;387.20056,633.89948;687.35120,829.62897;654.27692,823.89740;516.69287,738.75641;390.10764,621.99426;518.40881,719.58777;356.74579,671.57617;754.94531,845.78400;484.41171,724.67950;399.71805,702.71484;378.56238,637.62756;435.71350,628.17084;559.14673,53.482506;417.99997,725.85992;532.48706,711.29749;352.17267,672.54388;421.95151,628.86914;336.40613,668.55469;690.80817,852.99695;647.04260,847.00653;539.24628,119.09904;434.23328,738.67908;384.29330,62.580700;862.97699,266.69589;455.29346,47.657482;405.70880,653.27307;587.49677,112.30902;460.46634,796.10461;570.01923,123.93132;430.36334,653.94727;529.82129,74.714424;187.76447,540.92163;387.12000,31.590855;397.48755,54.040005;907.65363,170.78915;1059.4287,228.77196;327.05548,816.38916;603.25342,49.327290;863.31897,227.75424;348.98474,792.27228;857.96533,199.08241;1084.7850,233.81955;489.66577,595.14673];
 Point1 =[724.71436,850.77020;371.11609,664.90375;547.80957,764.37512;525.17426,772.77936;505.18951,761.67633;425.99838,665.37115;530.03485,726.89929;414.53741,851.45892;385.14536,660.84906;695.45490,840.09772;662.65839,836.46918;519.79675,758.75171;387.04755,648.69403;521.37524,739.22742;357.16089,700.21326;763.66351,852.91034;487.02643,746.51904;402.69095,729.35889;376.83527,665.56537;433.51724,651.96576;521.92346,68.013618;421.05270,750.98572;533.94574,730.09155;353.05063,701.50122;419.60675,654.31189;336.42163,699.05231;699.84955,863.16663;654.86865,859.70398;505.85379,134.82147;437.98547,763.06458;347.20193,87.611122;839.36523,264.87683;417.54456,68.360184;404.38141,679.30249;553.48370,125.27893;467.30707,819.36584;537.64807,138.13228;429.04501,678.50696;493.45453,91.087784;180.67621,578.02478;347.22858,56.859200;359.40833,78.548256;879.02948,166.10066;1034.3910,215.67024;336.67657,845.24286;565.31647,61.546837;837.62415,225.60124;355.89005,821.37708;831.12457,198.77525;1059.7500,219.32419;485.15564,616.58722];
 
%% Initialize: covariance matrix Ph observation equation
h = 1232; % size of image
w = 1640;
%-- Separate the location of points into 9 cells
counter = zeros(9,1);
m_o = cell(9,1); % orgianl match points in reference image
m_p = cell(9,1); % prime match points in alternative image
for i = 1:size(Point1,1)
    if Point1(i,1) <= fix(h/3) && Point1(i,2) <= fix(w/3) &&...
       Point2(i,1) <= fix(h/3) && Point2(i,2) <= fix(w/3) 
        counter(1) = counter(1)+1;
        m_o{1}(counter(1),:) = Point1(i,:);
        m_p{1}(counter(1),:) = Point2(i,:);
    elseif Point1(i,1) <= fix(2*h/3) && Point1(i,1)>fix(1*h/3) && Point1(i,2) <= fix(w/3) &&...
           Point2(i,1) <= fix(2*h/3) && Point2(i,1)>fix(1*h/3) && Point2(i,2) <= fix(w/3) 
        counter(2) = counter(2)+1;
        m_o{2}(counter(2),:) = Point1(i,:);
        m_p{2}(counter(2),:) = Point2(i,:);
    elseif Point1(i,1) <= fix(h) && Point1(i,1)>fix(2*h/3) && Point1(i,2) <= fix(w/3) &&...
           Point2(i,1) <= fix(h) && Point2(i,1)>fix(2*h/3) && Point2(i,2) <= fix(w/3)
        counter(3) = counter(3)+1;
        m_o{3}(counter(3),:) = Point1(i,:);
        m_p{3}(counter(3),:) = Point2(i,:);
    elseif Point1(i,1) <= fix(h/3) && Point1(i,1)>0 && Point1(i,2) <= fix(2*w/3) && Point1(i,2) > fix(w/3) &&...
           Point2(i,1) <= fix(h/3) && Point2(i,1)>0 && Point2(i,2) <= fix(2*w/3) && Point2(i,2) > fix(w/3) 
        counter(4) = counter(4)+1;
        m_o{4}(counter(4),:) = Point1(i,:);
        m_p{4}(counter(4),:) = Point2(i,:);
    elseif Point1(i,1) <= fix(2*h/3) && Point1(i,1)>fix(1*h/3) && Point1(i,2) <= fix(2*w/3) && Point1(i,2) > fix(w/3) &&...
           Point2(i,1) <= fix(2*h/3) && Point2(i,1)>fix(1*h/3) && Point2(i,2) <= fix(2*w/3) && Point2(i,2) > fix(w/3) 
        counter(5) = counter(5)+1;
        m_o{5}(counter(5),:) = Point1(i,:);
        m_p{5}(counter(5),:) = Point2(i,:);    
    elseif Point1(i,1) <= fix(h) && Point1(i,1)>fix(2*h/3) && Point1(i,2) <= fix(2*w/3) && Point1(i,2) > fix(w/3) &&...
           Point2(i,1) <= fix(h) && Point2(i,1)>fix(2*h/3) && Point2(i,2) <= fix(2*w/3) && Point2(i,2) > fix(w/3)  
        counter(6) = counter(6)+1;
        m_o{6}(counter(6),:) = Point1(i,:);
        m_p{6}(counter(6),:) = Point2(i,:);
    elseif Point1(i,1) <= fix(h/3) && Point1(i,1)>0 && Point1(i,2) <= fix(w) && Point1(i,2) > fix(2*w/3) &&...
           Point2(i,1) <= fix(h/3) && Point2(i,1)>0 && Point2(i,2) <= fix(w) && Point2(i,2) > fix(2*w/3)
        counter(7) = counter(7)+1;
        m_o{7}(counter(7),:) = Point1(i,:);
        m_p{7}(counter(7),:) = Point2(i,:);
    elseif Point1(i,1) <= fix(2*h/3) && Point1(i,1)>fix(1*h/3) && Point1(i,2) <= fix(w) && Point1(i,2) > fix(2*w/3) &&...
           Point2(i,1) <= fix(2*h/3) && Point2(i,1)>fix(1*h/3) && Point2(i,2) <= fix(w) && Point2(i,2) > fix(2*w/3) 
        counter(8) = counter(8)+1;
        m_o{8}(counter(8),:) = Point1(i,:);
        m_p{8}(counter(8),:) = Point2(i,:);    
    elseif Point1(i,1) <= fix(h) && Point1(i,1)>fix(2*h/3) && Point1(i,2) <= fix(w) && Point1(i,2) > fix(2*w/3) &&...
            Point2(i,1) <= fix(h) && Point2(i,1)>fix(2*h/3) && Point2(i,2) <= fix(w) && Point2(i,2) > fix(2*w/3)
        counter(9) = counter(9)+1;
        m_o{9}(counter(9),:) = Point1(i,:);
        m_p{9}(counter(9),:) = Point2(i,:);    
    end
end

ido = cellfun('length',m_o);
m_o(ido==0)=[];
idp = cellfun('length',m_p);
m_p(idp==0)=[];

%-- compute covariance matrix for each cell C_x and Jacobian J = dX/dh
len_m = size(m_o,1);
k = 1;
while k<=len_m
    if size(m_o{k},1) == 1
        m_o(k)=[];
        m_p(k)=[];
        len_m = len_m-1;
    else
        k = k+1;
    end
end

covariance_xi = cell(len_m,1);
covariance_x = [];
J_h = [];
error_x2 = [];
H0(3,1) = T2(3,1);
H0(3,2) = T2(3,2);
% H0(1,3) = T2(1,3)*1e-5;
% H0(2,3) = T2(2,3)*1e-5;
for i = 1:len_m
    
    p1 = ones(3,size(m_o{i},1));
    
    P1 = m_o{i}' ;
    P2 = m_p{i}' ;

    p1(1:2,:) = P1;
    p2 = H0*p1;
    
    P2_p = p2(1:2,:);
    errorm = P2-P2_p;
    errorv = reshape(errorm,[],1);
    error_x2 = [error_x2;errorv]; 
    
    covariance_xi{i} = cov(errorm');
    
    for l = 1:size(m_o{i},1)  
        covariance_x = blkdiag(covariance_x,covariance_xi{i});
    end
    
 %-- compute Jacobian J   
    m1 = p1;
    m2 = ones(3,size(m_o{i},1));
    m2(1:2,:) = P2_p;
    for j = 1:size(m_o{i},1)
        J_hj = [m1(1,j) m1(2,j) m1(3,j) 0 0 0 -m2(1,j)*m1(1,j)/p2(3,j) -m2(1,j)*m1(2,j)/p2(3,j) -m2(1,j)*m1(3,j)/p2(3,j);
               0 0 0 m1(1,j) m1(2,j) m1(3,j) -m2(2,j)*m1(1,j)/p2(3,j) -m2(2,j)*m1(2,j)/p2(3,j) -m2(2,j)*m1(3,j)/p2(3,j)]/p2(3,j);
        J_h = [J_h;J_hj];  
    end
end

h = reshape(H0',[9,1]);
v = h-sqrt(sum(h.*h))*[0 0 0 0 0 0 0 0 1]';

A = eye(9)-2*v*v'/(v'*v);

A = A(:,1:8);

S = A*inv(A'*J_h'*inv(covariance_x)*J_h*A)*A';
S = S/norm(S,2);
P = S;  

%-- Rearrange SURF points

Point1 = [];
Point2 = [];

for i = 1:len_m
    Point1 = [Point1;m_o{i}];
    Point2 = [Point2;m_p{i}];
end
global PP1 PP2
PP1 = Point1' ;
PP2 = Point2' ;
%% Define process equation and observation equation
Process = @(h)[h(1);h(2);h(3);h(4);h(5);h(6);h(7);h(8);h(9)];
Y = Observe(h);
[~,errm] = CovarianceEstPoints(Y);
% Y = Y+errm;
P = CovarianceHomo(h,covariance_x);
Q = 0.00*ones(9,9);
R = covariance_x;

%% UKF
 % [x,P]=ukf(Process,h,P,Observe,Y,Q,R); 
n=9; 
N=1;                                     % total dynamic steps
error_all = zeros(N,1);
for k=1:N
    [h,P]=ukf(Process,h,P,@Observe,reshape(PP2,[],1),Q,R);   % ukf     
    P = 0.5*P+0.5*P';
    P = P + 0.01*eye(9);
    h(3) = h(3)*1e-5;
    h(6) = h(6)*1e-5;
    Y = Observe(h);                     % measurments
    h = Process(h);                     % update process
    [R,errm] = CovarianceEstPoints(Y);
%     Y = Y+errm;
    errorm = reshape(PP2,[],1)-Y;
    error_all(k) = mean(abs(errorm));
    
end 
h = h/h(9); 
H = reshape(h,[3,3])'
error_all 




%% Compute Kalman gain & Estimate the homography
% h = reshape(H0,[9,1]);
% 
% error_all = mean(abs(error_x2));
% RMSE = 0*mean(abs(error_x2));
% 
% for k = 1:50
%     
% 
%     H = reshape(h,[3,3])';
%     error_x2 = [];
%     
%     J_h = [];
%     J_x = [];
%     Mdistance = [];
%     
%     p1 = ones(3,size(Point1,1));
%     
%     P1 = Point1' ;
%     P2 = Point2' ;
% 
%     p1(1:2,:) = P1;
%     p2 = H*p1;
% 
%     P2_p = p2(1:2,:);
%         
%         
%         
%     %-- compute Jacobian J   
%     m1 = p1;
%     m2 = ones(3,size(Point1,1));
%     m2(1:2,:) = P2_p;
%     for j = 1:size(P1,2)
%             
%         errorm = P2(:,j)-P2_p(:,j);
%         
%         
%         errorv = reshape(errorm,[],1);
%         
%         
%         error_x2 = [error_x2;errorv]; 
%         
%             
%         J_hj = [m1(1,j) m1(2,j) m1(3,j) 0 0 0 -m2(1,j)*m1(1,j)/p2(3,j) -m2(1,j)*m1(2,j)/p2(3,j) -m2(1,j)*m1(3,j)/p2(3,j);
%                 0 0 0 m1(1,j) m1(2,j) m1(3,j) -m2(2,j)*m1(1,j)/p2(3,j) -m2(2,j)*m1(2,j)/p2(3,j) -m2(2,j)*m1(3,j)/p2(3,j)]/p2(3,j);
%         J_h = [J_h;J_hj];  
%             
%             
%         J_xj = [h(1)-m2(1,j)*h(7)/p2(3,j) h(2)-m2(1,j)*h(8)/p2(3,j)  h(3)-m2(1,j)*h(9)/p2(3,j);
%                     h(4)-m2(2,j)*h(7)/p2(3,j) h(5)-m2(2,j)*h(8)/p2(3,j)  h(6)-m2(2,j)*h(9)/p2(3,j)]/p2(3,j);
%         J_x = [J_x;J_xj];
%             
%             
%         Q_j = J_xj*Lambda*J_xj' + J_hj*S*J_hj';
%         d_j = errorm'*pinv(Q_j)*errorm;
%             
%         Mdistance = [Mdistance;d_j];
%         
%             
%     end
%       
%     M = J_h;
%     W = J_x*Lambda*J_x';
%     L = W+M*S*M';
%     K = S* M'*pinv(W+M*S*M');
%     h = h + K*error_x2;
%     Ke = K*error_x2;
%     
%     S = (eye(9)-K * M)*S;
%     
%     S = S/norm(S,2);
%     
%     error_all(k) = mean(abs(error_x2));
%     P2_reshape = reshape(P2,[],1);
%     RMSE(k)= mean(error_x2.^2./P2_reshape);
%     
%     he = h/h(9);
%     H = reshape(he,[3,3]);
%     
%     %--Exclude outliers
%     threshold = (max(Mdistance)-min(Mdistance))*0.9+min(Mdistance);    
%     a = 1;
%     len_Md = size(Mdistance,1);
%     while k>=1 && a <= len_Md-1
%         if Mdistance(a)>threshold
%             Point1(a,:)  = [];
%             Point2(a,:)  = [];
%             Mdistance(a) = [];
%             len_Md = len_Md-1;
%         else
%             a = a+1;
%         end     
%     end
%     
% %     tform_bar = projective2d(H);
% %     
% %     [distorted_R,distorted_G1,distorted_G2,distorted_B] = Separate3Planes(distorted);
% %     cb_ref = imref2d(size(original_gs));
% %     recovered_R  = imwarp(distorted_R,tform_bar,'cubic','OutputView',cb_ref);
% %     recovered_G1 = imwarp(distorted_G1,tform_bar,'cubic','OutputView',cb_ref);
% %     recovered_G2 = imwarp(distorted_G2,tform_bar,'cubic','OutputView',cb_ref);
% %     recovered_B  = imwarp(distorted_B,tform_bar,'cubic','OutputView',cb_ref);
% %     recovered = Combine3Planes(recovered_R,recovered_G1,recovered_G2,recovered_B);
% %     
% %     figure
% %     imshowpair(original,recovered,'falsecolor')
%     
% end
% % figure 
% % 
% % plot(RMSE)
% % 
% % xlabel('iteration')
% % ylabel('relative mean squared error')
% 
% figure
% 
% plot(error_all)
% 
% xlabel('Iteration','FontSize',16)
% ylabel('Absolute error','FontSize',16)
% H
% error_all(end)
% 
% % H = T2;
% %H = T3;
% tform_bar = projective2d(H);
%     
% [distorted_R,distorted_G1,distorted_G2,distorted_B] = Separate3Planes(distorted);
% cb_ref = imref2d(size(original_gs));
%     recovered_R  = imwarp(distorted_R,tform_bar,'cubic','OutputView',cb_ref);
%     recovered_G1 = imwarp(distorted_G1,tform_bar,'cubic','OutputView',cb_ref);
%     recovered_G2 = imwarp(distorted_G2,tform_bar,'cubic','OutputView',cb_ref);
%     recovered_B  = imwarp(distorted_B,tform_bar,'cubic','OutputView',cb_ref);
%     recovered = Combine3Planes(recovered_R,recovered_G1,recovered_G2,recovered_B);
%     img_recovered(altNum) = {recovered*255};
%     if error_all(end)>=10
%         H_kalman(altNum) ={eye(3)};
%     else
%         H_kalman(altNum) = {H};
%     end
% %     figure
% %     imshowpair(original,recovered,'falsecolor')    
% %     
% end
% 
% len_H = size(H_kalman,2);
% k = 1;
% while k<=len_H
%     temp_H1 = H_kalman(k);
%     temp_H2 = cell2mat(temp_H1);
%     if temp_H2==eye(3)
%         img_recovered(k) = [];
%         H_kalman(k) = [];
%         len_H = len_H-1;
%     else
%         k = k+1;
%     end
% end
% N = size(img_recovered,2)+1;

%% MSE
% [m,n] = size(original);
%Z = original-recovered;
%MSE = 1/m*1/n*sum(sum(Z.*Z))
%PSNR = 10*log10(255^2/MSE)
% 
