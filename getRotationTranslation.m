function [angle,rot,trans] = getRotationTranslation(acc,gyro,timeStamp,imgStamp)
%load sensorData

load('noise_inertia.mat')

acc(:,1) = acc(:,1)-mean_acc(1);
acc(:,2) = acc(:,2)-mean_acc(2);
acc(:,3) = acc(:,3)-mean_acc(3);


gyro(:,1) = gyro(:,1)-mean_gyro(1);
gyro(:,2) = gyro(:,2)-mean_gyro(2);
gyro(:,3) = gyro(:,3)-mean_gyro(3);
imgNum = length(imgStamp);
rot = cell(imgNum,1); % current rotation matrix in the initial frame
trans = cell(imgNum,1);
angle = cell(imgNum,1);
acc = acc*0;
gyro = gyro*1+0.0000;
delay = 0;
imgStamp = imgStamp*1;

for n = 1:imgNum-1
    R = eye(3);
    v = zeros(3,1);
    T = zeros(3,1); 
    if imgStamp(n+1)-imgStamp(n)>10
        imgStampv=imgStamp(n)+8;
    else
        imgStampv=imgStamp(n+1);
    end

    for i = imgStamp(n)+1+delay:imgStampv
        R0 = R;
        v0 = v;
        T0 = T;
        
        dt = (timeStamp(i+1)-timeStamp(i))/1e9; %unit second
        dphi = gyro(i,:)*dt;
        
        R = [1 -dphi(3) dphi(2); dphi(3) 1 -dphi(1); dphi(2) dphi(1) 1]*R0;
        
        if i == imgStamp(n)+1+delay 
            v = v0 + inv(R)*acc(i,:)'*dt/2;
        else
            v = v0 + (inv(R0)*acc(i-1,:)' + inv(R)*acc(i,:)'*dt/2);
        end
        
        T = T0 + (v0 + v)*dt/2;
    end
    eulerangle = convertR2E(R);
    T;
    rot{n} = R;
    trans{n} = T;
    angle{n} = convertR2E(R);
end

end

function eulerangle = convertR2E(R)
   eulerangle = zeros(3,1);
   eulerangle(1,1) = atan(R(3,2)/R(3,3)); % euler angle for x-axis
   eulerangle(2,1) = asin(-R(3,1)); % euler angle for y-axis
   eulerangle(3,1) = atan(R(2,1)/R(1,1)); % euler angle for z-axis
end
