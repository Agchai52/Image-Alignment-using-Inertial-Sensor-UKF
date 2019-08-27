%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Main_Lowlight                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
%% Initialize
N = 14;
img_offset  = cell(1,N);
img_aligned = cell(1,N);


%% Simulation
% img0 = imread('Peper_reference.png');
% img0gs = rgb2gray(img0);
% sigma = 0.1;
% I0 = imnoise(I0/255,'gaussian',0,sigma^2)*255;
% for i = InitialFrame+1:N-1
%     I1 = PrepareOneOffsetImage(I0,sigma);
%     img_offset(i) = {I1};
% end



%% Real Images

FileName = '20180125_154100_0.dng';%'20171101_172811_0.dng';%'20171113_160710_0.dng';'20180124_145422_0.dng''20171105_201630_0.dng'
FileDate = FileName(1:16);
Format = '.dng';

% Read Images
I0 = double(imread(FileName)); % as the reference

for i = 1:N-1
    indexnum = num2str(i);
    FileName = strcat(FileDate,indexnum,Format);
    img_offset(i) = {double(imread(FileName))};
end
% img_aligned=img_offset;
%% Alignment by Inertia Sensor Data
% % % % img_aligned = img_recovered; 
run Alignment_InertiaSensor

img_offset = img_recovered;
%img_aligned=img_offset;
validFrameNum = size(img_recovered,2)-InitialFrame+1
%% Alignment by Gaussian Pyramid\
run Alignment_GaussianPyramid

%% Merge 
% Seperate Planes
[R0,G10,G20,B0] = Separate3Planes(I0);
R = cell(1,N);
G1= cell(1,N);
G2= cell(1,N);
B = cell(1,N);
for i = 1:N-1
    [tempR,tempG1,tempG2,tempB] = Separate3Planes(cell2mat(img_aligned(i)));
    R(i)  = {tempR};
    G1(i) = {tempG1};
    G2(i) = {tempG2};
    B(i)  = {tempB};
end

% Merge each plane
disp('Merging Plane R');
R_merged = MergeAllImages(R0,R);
disp('Merging Plane G1');
G1_merged = MergeAllImages(G10,G1);
disp('Merging Plane G2');
G2_merged = MergeAllImages(G20,G2);
disp('Merging Plane B');
B_merged = MergeAllImages(B0,B);

% Combine planes
img_merged = Combine3Planes(R_merged,G1_merged,G2_merged,B_merged);
disp('Merge Ends');
%% Comparison
%I_wn = wiener2(I0,[5,5]);
%figure,imshow(I0/255)
figure,imshow(img_merged/255)
%figure,imshow(I_wn/255)

Photo=demosaic(im2uint16(img_merged/255),'rggb');
figure,imshow(Photo,[])

run ProcessMergedImage
% ssimval_Noise = ssim(I0,I00)
% ssimval_Denoise = ssim(img_merged,I00)
% 
% 
% psnr_Noise = mpsnr(I0,I00)
% psnr_Denoise = mpsnr(img_merged,I00)