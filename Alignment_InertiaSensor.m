clc
clear all
close all

% #Image
N = 14;
InitialFrame = 0;
FileName = '20171101_172811_0.dng';%'20171101_172811_0.dng';%'20171105_201630_0.dng''20171105_171221_0.dng''20180125_154100_0.dng'
FileDate = FileName(1:16);
Format = '.dng';
run SURF_Homography
run ReadSensorData


img_recovered = cell(1,N-1);
H_kalman = cell(1,N-1);
for altNum = 1:N-1
    T2 = H_SURF{altNum} %image Current->Target
    sign_T2 = sign(reshape(T2',[9,1]));
    Point2 = SURF_Point2{altNum};
    Point1 = SURF_Point1{altNum};   

    T2i = K_bar*inv(T2)*inv(K_bar); %camera Target-> Current
    T2i = T2i/T2i(3,3)

    % RelativeCameraPose
    angle_bar = angle{altNum};
    R_bar = rot{altNum};
    t_bar = trans{altNum};
    H_bar = R_bar%-[0 0 0; 0 0 0;t_bar']; %camera Target-> Current
    H_bar = R_bar*(eye(3)+[0;0;1]*t_bar');
    
    H0 = inv(K_bar)*inv(H_bar)*(K_bar);  %image Current->Target

    H0 = H0/H0(3,3);
    T3 = H0;
    
%     H0(3,1) = T2(3,1);
%     H0(3,2) = T2(3,2);
%     H0(1,3) = T2(1,3);
%     H0(2,3) = T2(2,3);
    
H0;


% H0 = T2;
% Point1 = [1016.6059,453.36786;1266.1906,884.72424;1417.8013,965.51819;808.56848,440.00635;707.53571,413.95816;1120.1685,709.35565;478.26285,175.65819;908.69885,391.47040;321.26364,88.399353;749.43243,331.84995;625.03882,43.610035;1262.7933,890.78894;886.63129,104.45791;1033.7781,453.11569;422.25003,399.66568;263.44461,56.629620;449.82339,96.700310;816.80896,625.68018;624.16541,428.53552;927.89398,61.459988];
% Point2 = [1016.7157,452.62393;1266.3505,884.43237;1417.9493,964.75031;808.65472,439.24149;707.93915,413.08551;1120.0406,708.58069;478.91443,174.89153;908.15735,390.15512;321.65485,87.655716;749.60980,330.66812;625.14636,42.380062;1262.7843,890.15228;886.88812,103.74699;1033.8928,452.64990;422.47192,398.73096;264.01556,55.404697;450.33740,95.671310;816.95905,625.02155;624.59473,427.36212;928.12469,60.582428];
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
    

end

h = reshape(H0',[9,1]);

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
[covariance_x,errm] = CovarianceEstPoints(Y);
% Y = Y+errm;
P = CovarianceHomo(h,covariance_x);
Q = 0.000*ones(9,9);
R = covariance_x;

%% UKF
 % [x,P]=ukf(Process,h,P,Observe,Y,Q,R); 
n=9; 
StepN=1;                                     % total dynamic steps
error_all = zeros(StepN,1);
for k=1:StepN
    [h,P]=ukf(Process,h,P,@Observe,reshape(PP2,[],1),Q,R);   % ukf     
    P = 0.5*P+0.5*P';
    P = P + 0.01*eye(9);
    h = sqrt(h.*conj(h));
    h(3) = h(3)*1e-5;
    h(6) = h(6)*1e-5;
    sign_h = sign(h);
    h = sign_T2.*sign_h.*h;
    H = H0;
    h = reshape(H',[9,1]);
    Y = Observe(h);                     % measurments
    h = Process(h);                     % update process                        
%     h = reshape(H',[9,1]);
    [R,errm] = CovarianceEstPoints(Y);
%     Y = Y+errm;
    
    
    errorm = reshape(PP2,[],1)-Y;
    error_all(k) = mean(abs(errorm));
    
end 
h = h/h(9); 
H = reshape(h,[3,3])'

error_all(end) 
%% Compute Kalman gain & Estimate the homography
% figure 
% 
% plot(RMSE)
% 
% xlabel('iteration')
% ylabel('relative mean squared error')

% figure
% 
% plot(error_all)
% 
% xlabel('Iteration','FontSize',16)
% ylabel('Absolute error','FontSize',16)



%  H = T2; H = T3;
% H = H0;
tform_bar = projective2d(H);
distorted = distorted_all{altNum};   
[distorted_R,distorted_G1,distorted_G2,distorted_B] = Separate3Planes(distorted);
cb_ref = imref2d(size(original_gs));
    recovered_R  = imwarp(distorted_R,tform_bar,'cubic','OutputView',cb_ref);
    recovered_G1 = imwarp(distorted_G1,tform_bar,'cubic','OutputView',cb_ref);
    recovered_G2 = imwarp(distorted_G2,tform_bar,'cubic','OutputView',cb_ref);
    recovered_B  = imwarp(distorted_B,tform_bar,'cubic','OutputView',cb_ref);
    recovered = Combine3Planes(recovered_R,recovered_G1,recovered_G2,recovered_B);
    img_recovered(altNum) = {recovered*255};
    if error_all(end)>=70
        H_kalman(altNum) ={eye(3)};
    else
        H_kalman(altNum) = {H};
    end
    figure
    imshowpair(original,recovered,'falsecolor')    
    
end

len_H = size(H_kalman,2);
k = 1;
while k<=len_H
    temp_H1 = H_kalman(k);
    temp_H2 = cell2mat(temp_H1);
    if temp_H2==eye(3) 
        img_recovered(k) = [];
        H_kalman(k) = [];
        len_H = len_H-1;
    else
        k = k+1;
    end
end
N = size(img_recovered,2)+1;


