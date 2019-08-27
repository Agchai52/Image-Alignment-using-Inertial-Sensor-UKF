function [var] = VarDetect(T0,T)

T = T-T0;

T2 = T.*T;

RMS2 = mean2(T2);
Mean2 = mean2(T)^2;

var = (RMS2-Mean2);

end

