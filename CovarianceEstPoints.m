function [covariance_x,errm] = CovarianceEstPoints(Y)

global PP2
PointsNum = size(PP2,2);

P2_p = reshape(Y,[2,PointsNum]);
P2 = PP2;
        
errorm = P2-P2_p;
% errorm = Y;
errm_xi = mean(errorm')';
covariance_xi = cov(errorm'); 
covariance_x = [];
errm = [];
for l = 1:PointsNum  
    covariance_x = blkdiag(covariance_x,covariance_xi);
    errm = [errm;errm_xi];
end        

end