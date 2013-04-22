clc
clear all
close all

% *************************************************************************
% Heeger (1988).
%   Estimate optic flow for a rotating random texture and display the
%   estiamted optic flow together with the likelihood values for motions.
%
%   Copyright (C) 2013  Florian Raudies, 01/14/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.
% *************************************************************************

% Generate stimulus.
addpath('./../'); % add root directory to path to access auxiliary functions.
opt.h       = 128; % set parameters for stimulus.
opt.w       = 128;
opt.fNum    = 15;
% opt.omega   = pi/300;
% [ImgSeq maxSpeed] = randomTextureRotation(opt);
opt.expRate = 1.005;
[ImgSeq maxSpeed] = randomTextureExpansion(opt);
rmpath('./../');

% Set parameters for estimation method.
opt.sigmaV = 5*10^-1;
[Dx Dy L] = estimateOpticFlow2D(ImgSeq,opt);

% Display the estimated optic flow.
h       = size(ImgSeq,1);
w       = size(ImgSeq,2);
[Y X]   = ndgrid(1:h, 1:w); % pixel coordinates.
sample  = 8;
IndexX  = 1:sample:w;
IndexY  = 1:sample:h; 
len     = sample/maxSpeed;
figure('Position',[50 50 600 600]); 
quiver(X(IndexY,IndexX),      Y(IndexY,IndexX),...
       Dx(IndexY,IndexX)*len, Dy(IndexY,IndexX)*len,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));

% Display the motion likelihood values by construction of one image that
% contains likelihoods for all spatial positions as a tile. 
vyNum = size(L,3);
vxNum = size(L,4);
Img = ones([h*vyNum+(vyNum-1), w*vxNum+(vxNum-1)]);
for ivy = 1:vyNum,
    for ivx = 1:vxNum,
        Img((ivy-1)*(h+1)+(1:h),(ivx-1)*(w+1)+(1:w)) = squeeze(L(:,:,ivy,ivx));
    end
end
figure('Position',[650 50 600 600]);
% Low-pass filter to avoid alaising that lines show in 'imshow'.
Img = imfilter(Img,fspecial('gaussian',[5 5],1),'same','replicate');
imshow(Img,[0 1]);
title(['Each tile represents the motion likelihood for one velocity and'...
      'all spatial locations.']);
