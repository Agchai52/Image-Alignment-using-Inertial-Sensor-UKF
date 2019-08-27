%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Read sensor data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc
% clear all
% close all
% FileDate = '20171107_162543_';% steady error files '20171107_162543_'
facc = fopen([FileDate,'acc.txt'],'r');
acc = fscanf(facc,'%f %f %f %f\n',[4,Inf]);
fclose(facc);
acc = acc';

fgyro = fopen([FileDate,'gyro.txt'],'r');
gyro = fscanf(fgyro,'%f %f %f %f\n',[4,Inf]);
fclose(fgyro);
gyro = gyro';

fparam = fopen([FileDate,'param.txt'],'r');
param = fscanf(fparam,'%f %f %f %f\n');
fclose(fparam);
param = param';

timeStamp = acc(:,1);
acc = acc(:,2:4);
gyro = gyro(:,2:4);

acc_temp = acc;
acc(:,1) = acc_temp(:,2);
acc(:,2) = acc_temp(:,1);

gyro_temp = gyro;
gyro(:,1) = gyro_temp(:,2);
gyro(:,2) = gyro_temp(:,1);
gyro(:,3) = -gyro(:,3);

% mean_acc  = [mean(acc(:,1)),mean(acc(:,2)),mean(acc(:,3))]'
% mean_gyro = [mean(gyro(:,1)),mean(gyro(:,2)),mean(gyro(:,3))]'
% 
% cov_acc  = diag([cov(acc(:,1)),cov(acc(:,2)),cov(acc(:,3))])
% cov_gyro = diag([cov(gyro(:,1)),cov(gyro(:,2)),cov(gyro(:,3))])
% 
% save noise_inertia mean_acc mean_gyro cov_acc cov_gyro
% load noise_inertia




imgNum = N;%3
imgStamp = param(1:2*imgNum);
imgStamp = fix(downsample(imgStamp,2));
exposure = param(2*(imgNum+1)+3);
focal = param((imgNum+1)*2+1);
f = 3272/4*focal/4.536; % focal length in pixels focal_pixel = (focal_mm / sensor_width_mm) * image_width_in_pixels
% camera intrinsic matrix
h = 2464;
w = 3280;
% [h,w,~] = size(Image);
cp = [floor(w/2),floor(h/2)];
K_bar = [f 0 cp(1);0 f cp(2); 0 0 1]';
%  K_bar = [f 0 0;0 f 0; 0 0 1];

save sensorData timeStamp acc gyro imgStamp K_bar
% [angle,rot,trans] = getRotationTranslation(acc,gyro,timeStamp,imgStamp);
% run Example03Inertial
[rot,trans,angle] = InertiaSensorDataOutput(imgNum,H_SURF);


% angle_bar = angle{altNum};
% R_bar = rot{altNum};
% t_bar = trans{altNum};
% H_bar = R_bar - [0 0 0; 0 0 0;-t_bar']';

