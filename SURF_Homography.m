%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             SURF_Homography                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To run this, run Alignment_InertiaSensor
H_SURF = cell(1,N-1);
SURF_Point1 = cell(1,N-1);
SURF_Point2 = cell(1,N-1);
distorted_all = cell(1,N-1);

% Read Images
original = double(imread(FileName))/255; % as the reference

% Sperate each image into 3 images
[original_R,original_G1,original_G2,original_B] = Separate3Planes(original);
original_RGB(:,:,1) = original_R;
original_RGB(:,:,2) = (original_G1 + original_G2)/2;
original_RGB(:,:,3) = original_B;
original_gs = rgb2gray(original_RGB);


%% Matching SURF features
% Detecting the matching SURF features
ptsOriginal  = detectSURFFeatures(original_gs,'MetricThreshold',10,'NumScaleLevels',3,'NumOctaves',10);%(original_gs,'NumScaleLevels',3,'NumOctaves',3);
% Extract features
[featuresOriginal,validPtsOriginal] = ...
            extractFeatures(original_gs,ptsOriginal,'SURFSize',128);
        

for altNum =1:N-1
    % Sperate each image into 3 images
    indexnum = num2str(altNum);
    FileName = strcat(FileDate,indexnum,Format);
    distorted = double(imread(FileName))/255; % as the alternatives
    %     figure
%     imshowpair(original,distorted,'falsecolor')
    
    [distorted_R,distorted_G1,distorted_G2,distorted_B] = Separate3Planes(distorted);
    distorted_RGB(:,:,1) = distorted_R;
    distorted_RGB(:,:,2) = (distorted_G1 + distorted_G2)/2;
    distorted_RGB(:,:,3) = distorted_B;
    distorted_gs = rgb2gray(distorted_RGB);
    
    % Detecting the matching SURF features
    ptsDistorted = detectSURFFeatures(distorted_gs,'MetricThreshold',10,'NumScaleLevels',3,'NumOctaves',10);%(distorted_gs,'NumScaleLevels',3,'NumOctaves',3);
    % Extract features
    [featuresDistorted,validPtsDistorted] = ...
            extractFeatures(distorted_gs,ptsDistorted,'SURFSize',128);
        
        
% figure
% imshow(original_gs);
% hold on;
% strongestPoints = validPtsOriginal.selectStrongest(10);
% strongestPoints.plot('showOrientation',true);
% 
% figure
% imshow(distorted_gs);
% hold on;
% strongestPoints = validPtsDistorted.selectStrongest(10);
% strongestPoints.plot('showOrientation',true);


    % Find candidate matches.
    indexPairs = matchFeatures(featuresOriginal,featuresDistorted,'MatchThreshold',1.0);
    % Find point locations from both images.
    matchedOriginal  = ptsOriginal(indexPairs(:,1));
    matchedDistorted = ptsDistorted(indexPairs(:,2));
    % Coordinate of points
    Point2 = matchedOriginal.Location;
    Point1 = matchedDistorted.Location;

    
    % Analyze the feature locations.
    [tform, inlierDistorted,inlierOriginal] = ...
            estimateGeometricTransform(matchedDistorted.Location,...
                matchedOriginal.Location,'projective','MaxNumTrials',3000);

    T2 = tform.T;

    H_SURF{altNum} = T2;
    SURF_Point1{altNum} = inlierDistorted;
    SURF_Point2{altNum} = inlierOriginal;
    distorted_all{altNum} = distorted;
end   

clear altNum