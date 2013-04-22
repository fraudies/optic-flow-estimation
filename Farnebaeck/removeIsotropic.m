function T = removeIsotropic(T)
% removeIsotropic
%
% T is assumed to be an MxNx3x3 tensor field. 
% For details, see Appendix G of Gunnar Farnebäck's thesis "Polynomial
% Expansion for Orientation and Motion Estimation".

T(:,:,[1 5 9]) = T(:,:,[1 5 9]) - repmat(sum(T(:,:,[1 5 9]), 3)/3, [1 1 3]);
p = T(:,:,1,1).*T(:,:,2,2) + T(:,:,1,1).*T(:,:,3,3) + ...
    T(:,:,2,2).*T(:,:,3,3) - T(:,:,1,2).^2 - T(:,:,1,3).^2 - T(:,:,2,3).^2;
q = T(:,:,1,1).*T(:,:,2,3).^2 + T(:,:,2,2).*T(:,:,1,3).^2 + ...
    T(:,:,3,3).*T(:,:,1,2).^2 - 2*T(:,:,1,2).*T(:,:,1,3).*T(:,:,2,3) - ...
    T(:,:,1,1).*T(:,:,2,2).*T(:,:,3,3);
beta        = sqrt(-4*p/3);
alpha       = real(acos(3*q./(-eps+p.*beta))/3);
lambda3     = beta.*cos(alpha+2*pi/3);
T(:,:,[1 5 9]) = T(:,:,[1 5 9]) - repmat(lambda3, [1 1 3]);
