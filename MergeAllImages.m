function img_merged = MergeAllImages(I0,I)
%% Robust Pairwise Temporal Merge
%N = 9; % number of reference and alternative image
c = 1e5; 
t = 32; % size of tile even
M = 13; % coef of spatial denoise filter
 
N = size(I,2);

[m0,n0] = size(I0);

padm = mod(m0,t);
padn = mod(n0,t);

m = m0 + padm;
n = n0 + padn;

InitialFrame = 0;

I0 = padarray(I0,[padm,padn],'replicate','post');
for i = InitialFrame+1:N-1
    tempI1 = cell2mat(I(i));
    tempI2 = padarray(tempI1,[padm,padn],'replicate','post');
    I(i) = {tempI2};
end

I0_tilt = zeros(m,n);
I_sd = zeros(m,n);
I_ot = zeros(m,n);

% Windows     
     % Analysis Window to Tf0_tilt
     n1 = -(t-0.5):(t-0.5);
     n2 = -(t-0.5):(t-0.5);
     [n1,n2] = meshgrid(n1,n2);
     window_an = cos(n1*pi/2/t).* cos(n2*pi/2/t);
     
     % Output Window
     n1 = 0:(2*t-1);
     n2 = 0:(2*t-1);
     [n1,n2] = meshgrid(n1,n2);
     window_ot = (0.5-0.5*cos(2*pi*(n1+0.5)/2/t)).*(0.5-0.5*cos(2*pi*(n2+0.5)/2/t));
      %cos(n1*pi/2/t)^2.*cos(n2*pi/2/t)^2;
      
% Coef of spatial denoise filter    
fsigma = 0:M/(2*t):(M-M/(2*t));
f = fsigma'*fsigma;
               
for i = 1:(t/2):(m-t) %129
    for j = 1:(t/2):(n-t) %89
        
        T0 = I0(i:i+t-1,j:j+t-1);  
        Tf0 = fft2(T0,2*t,2*t);
        Tf0_tilt = 1/N*Tf0;
        
        for index = 1:N-1
            tempI1 = cell2mat(I(index));
            
            T   = tempI1(i:i+t-1,j:j+t-1);
            Tf = fft2(T,2*t,2*t);
            var = VarDetect(T0,T);
            Af = abs(Tf0-Tf).^2./(abs(Tf0-Tf).^2+c*var);
            
            % Merge Equation 6
            % Tf0_tilt = 1/N*sum(Tfi+Afi.*(Tf0-Tfi)) i = 0...N
            Tf0_tilt = Tf0_tilt + 1/N*(Tf+Af.*(Tf0-Tf));
            
        end
       
 
      T0_tilt = real(ifft2(Tf0_tilt));
      T0_tilt = T0_tilt(1:t,1:t);
      I0_tilt(i:i+t-1,j:j+t-1) = T0_tilt(1:t,1:t);
     
    
     T11 = T0_tilt;
     T12 = flip(T11,2);
     T21 = flip(T11,1);
     T22 = flip(T21,2);
     
     T_rep=[T11,T12;T21,T22];
     
     Tf_rep = fft2(T_rep);
     Tf_an = fftshift(Tf_rep).* window_an;
      
     % Overlapped Titles
     Tf_ot = Tf_an.* window_ot;
     Tf_ot = ifftshift(Tf_ot);
     T_ot = real(ifft2(Tf_ot));
     Tf_ot = fft2(T_ot,2*t,2*t);
     I_ot(i:i+t-1,j:j+t-1) = T_ot(1:t,1:t);
    
   


     % Spatial Denoising
     
     varz = VarDetect(T0,T_ot(1:t,1:t));
     D_z  = Tf0-Tf_ot;
     A_z  = abs(D_z).^2./(abs(D_z).^2 + f*varz*8/N +1);
     
     Tf_sd = Tf_ot.*A_z;
   
     T_sd = real(ifft2(Tf_sd));
     I_sd(i:i+t-1,j:j+t-1) = T_sd(1:t,1:t);
     
    end
end


img_merged = I_sd(1:m0,1:n0);
end

