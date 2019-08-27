function Y = Observe(h)
global PP1 PP2
H = reshape(h,[3,3])';

p1 = ones(3,size(PP1,2));
    
p1(1:2,:) = PP1;
p2 = H*p1;

Y = p2(1:2,:);
% Y = Y-PP2;
Y = reshape(Y,[],1);

end

