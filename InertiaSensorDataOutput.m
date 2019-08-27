function [rot,trans,angle] = InertiaSensorDataOutput(imgNum,H_SURF)

load('noise_inertia.mat')
load('sensorData.mat')
acc(:,1) = acc(:,1)-mean_acc(1);
acc(:,2) = acc(:,2)-mean_acc(2);
acc(:,3) = acc(:,3)-mean_acc(3);


gyro(:,1) = gyro(:,1)-mean_gyro(1);
gyro(:,2) = gyro(:,2)-mean_gyro(2);
gyro(:,3) = gyro(:,3)-mean_gyro(3);

acc  = 1*acc;
gyro = 1*gyro;

% Xaug = [wx wy wz anglex angley anglez ax ay az vx vy vz tx ty tz]]
% n=15;      %number of state
q=0.00001;    %std of process 
r=0.0001;    %std of measurement
Qi = q^2*eye(3); 
Q  = blkdiag(zeros(3,3),Qi,zeros(3,3),Qi,Qi); % covariance of process
R  = blkdiag(q^2*eye(3),r^2*eye(3));        % covariance of measurement  

X0 = zeros(15,1);                               % initial state
Xaug_test = zeros(15,2); 
P_test = zeros(2,1);
P_t = zeros(3,2);
Pw = cov_gyro;
Pa = cov_acc;
P = blkdiag(Pw,zeros(3,3),Pa,zeros(6,6)); % initial state covraiance
P = 0.5*P+0.5*P';
P = P + 0.01*eye(15);
% imgNum = length(imgStamp);
rot = cell(imgNum,1); % current rotation matrix in the initial frame
trans = cell(imgNum,1);
angle = cell(imgNum,1);
imgStamp = imgStamp*1;

global dt
for n = 1:imgNum-1
    Xaug(1:9,1) = zeros(9,1); 
    Xaug(13:15,1) = zeros(3,1); % only keep v
    if imgStamp(n+1)-imgStamp(n)>10
        imgStampv=imgStamp(n)+10;
    else
        imgStampv=imgStamp(n+1);
    end
    H_SURFInv = H_SURF{n};
    H_SURFRea = K_bar*inv(H_SURFInv)*inv(K_bar);
    H_SURFRea = H_SURFRea/H_SURFRea(3,3)
    [R1,T1,t1,angle1,R2,T2,t2,angle2] = decompHomoMatrix2(H_SURFRea);
%     lambda_t = [lambda1,lambda2];
    R_t{1} = R1;
    R_t{2} = R2;
    angle_t{1} = angle1;
    angle_t{2} = angle2;
    P_t(:,1) = t1; 
    P_t(:,2) = t2; 
    P_test = abs([R1(3,1)/R1(3,3)-H_SURFRea(3,1),R2(3,1)/R2(3,3)-H_SURFRea(3,1)]);
    node = find(P_test==min(P_test))
%     R_tt = R_t{node}
    t_t = P_t(:,node)
    z_t = [angle_t{node};P_t(:,node)];
    for i = imgStamp(1)+1:imgStamp(n+1)+1
        dt = (timeStamp(i+1)-timeStamp(i))/1e9; %unit second
        Xaug(1:3,1) = gyro(i,:);
        Xaug(7:9,1) = acc(i,:);
        z = (i-imgStamp(1))/(imgStamp(n+1)-imgStamp(1)+1)*z_t;
        
        [Xaug, P] = ukf(@f_process,Xaug,P,@h_measure,z,Q,R);
  
        P = 0.5*P+0.5*P';
        P = P + 0.01*eye(15);
        Xaug = f_process(Xaug);
    end
    eulerangle = Xaug(4:6,1);   
    T_bar = Xaug(13:15,1)%t_t%
    R_bar = convertE2R(eulerangle);
    H_bar = R_bar*(eye(3)+[0;0;1]*T_bar')
    rot{n} = R_bar;
    trans{n} = T_bar;
    angle{n} = eulerangle;
end