%=========================================================================%
%                          Process Merged Image                           %
%=========================================================================%
%% Parameter of LGE Nexus 5
% Raw image to be processed
% raw = double(imread('20180124_145422_0.dng'));%double(imread('Reference.tiff'));
raw = img_merged;
% For linearization
darkness = 64; % from decraw
saturation = 1023;

% For white balance
r_multiplier = 1.910774; % no white balance of camera 1.6;
b_multiplier = 1.758186; % no white balance of camera 2.8;
% r_multiplier = 1.242718; % use white balance of camera
% b_multiplier = 3.555556; % use white balance of camera

% For color correction
ColorMatrix1 = [0.6171875,-0.140625,-0.0390625;-0.4609375,1.1953125,0.2265625;-0.0859375,0.1953125,0.546875];
ColorMatrix2 = [1.3046875,-0.5625,0.0078125;-0.171875,1.0625,0.359375;0.0390625,0.0546875,0.6171875];
CameraCalibration1 = [0.984375,0,0;0,1,0;0,0,0.9609375];
CameraCalibration2 = [0.984375,0,0;0,1,0;0,0,0.9609375];
AsShotNeutral = [0.562500000000000 1 0.453125000000000];
CalibrationIlluminant1 = 21;
CalibrationIlluminant2 = 17;
ForwardMatrix1 = [0.7578125,0.25,-0.046875;0.296875,0.9765625,-0.2734375;0.0078125,-0.2265625,1.046875];
ForwardMatrix2 = [0.609375,0.2734375,0.078125;0.0625,1.09375,-0.15625;-0.109375,-0.3515625,1.28125];
W1 = 0.8; % picked myself
W2 = 1-W1;



%% Linearization
% darkness = min(raw(:));
% saturation = max(raw(:));
lin_bayer = (raw-darkness)/(saturation-darkness); % normalized into [0,1];
lin_bayer = max(0,min(lin_bayer,1)); % ensure no element outside [0,1];
%figure, imshow(lin_bayer)
%% White Balancing
wb_multipliers = [r_multiplier,1,b_multiplier]; % for particular condition, from dcraw;
mask = wbmask(size(lin_bayer,1),size(lin_bayer,2),wb_multipliers);
balanced_bayer = lin_bayer .* mask;
%figure, imshow(balanced_bayer)
% %% Sharpening
% h = [-1 1 -1;1 8 -1;-1 1 -1];
% nl_sharp = imfilter(balanced_bayer,h);
% figure, imshow(nl_sharp)
% balanced_bayer = nl_sharp;
%% Demosaicking
temp = uint16(balanced_bayer/max(balanced_bayer(:)) * (2^16-1));%(2^16-1)
lin_rgb = double(demosaic(temp,'rggb'))/(2^16-1);
%figure, imshow(lin_rgb)

%% CMOSCrop Info.DefaultCropSize [3272,2456] origin[8,8]
lin_crop=lin_rgb(0+8+1:size(lin_bayer,1)-8,0+8+1:size(lin_bayer,2)-8,:);

%% Color Space Conversion
AB = [1,0,0;0,1,0;0,0,1];
CM = W1*ColorMatrix1+W2*ColorMatrix2;
CC = W1*CameraCalibration1+W2*CameraCalibration2;
FM = W1*ForwardMatrix1+W2*ForwardMatrix2;
XYZ2Cam = AB*CC*CM;
% sRGB2XYZ is an unchanged standard
sRGB2XYZ = [0.4124564 0.3575761 0.1804375;0.2126729 0.7151522 0.0721750;0.0193339 0.1191920 0.9503041];
% Here XYZ2Cam is only for Nikon D3X, can be found in adobe_coeff in dcraw.c
%XYZ2Cam = [7171 -1986 -648;-8085 15555 2718;-2170 2512 7457]/10000;

sRGB2Cam = XYZ2Cam * sRGB2XYZ;
sRGB2Cam = sRGB2Cam./ repmat(sum(sRGB2Cam,2),1,3); % normalize each rows of sRGB2Cam to 1
Cam2sRGB = (sRGB2Cam)^-1;
lin_srgb = apply_cmatrix(lin_crop, Cam2sRGB);
lin_srgb = max(0,min(lin_srgb,1)); % Always keep image clipped b/w 0-1
%figure, imshow(lin_srgb)
%% Brightness and Gamma Correction
grayim = rgb2gray(lin_srgb); % Consider only gray channel
grayscale = 0.25/mean(grayim(:));
bright_srgb = min(1,lin_srgb * grayscale); % Always keep image value less than 1

nl_srgb = bright_srgb.^(1/1.95);

sharpen_srgb = imsharpen(nl_srgb,'Radius',2,'Amount',1.2);

[msp,nsp]=size(sharpen_srgb);
tic
nlevel = NoiseLevel(sharpen_srgb(1:fix(msp/2),1:fix(nsp/2)),7,0,0.9,1);
noiselevel = norm(mean(nlevel))*255
PSNR = 10*log10(255^2/noiselevel^2)
toc
% figure, imshow(nl_srgb(1200:1500,1800:2100,:))%(1700:2000,2700:3000,:))
figure, imshow(sharpen_srgb);%(1200:1500,1800:2100,:))%(1100:1400,1250:1550,:))(1000:1300,1650:1950,:)
% figure, imshow(nl_srgb)
%  figure, imshow(sharpen_srgb);
figure, imshow(sharpen_srgb(1200:1500,1800:2100,:))

%imwrite(nl_srgb,'img2.jpg','jpg')
%figure, imshow(nl_srgb(2000:2300,850:1200,:))
