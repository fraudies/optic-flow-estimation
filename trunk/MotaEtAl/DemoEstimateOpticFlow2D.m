clc
clear all
close all

% *************************************************************************
% Mota et al. (2001).
% It takes a few moments (~10sec) before the first flow is estimated and 
% displayed.
%   Copyright (C) 2013  Florian Raudies, 01/08/2013, Boston University.
%   License, GNU GPL, free software, without any warranty.
% *************************************************************************

% Generate the first two image sequences.
addpath('./../'); % add root directory to path to access auxiliary functions.
opt.h       = 128; % Set parameters for stimulus.
opt.w       = 128;
opt.fNum    = 11;
opt.Vel1    = [1 1];
opt.Vel2    = [0 -1];
[ImgSeq1 maxSpeed1] = randomTextures2linearMotions(opt);
opt.expRate = 1.02;
[ImgSeq2 maxSpeed2] = randomTexturesExpCon(opt);
opt.omega   = pi/100;
[ImgSeq3 maxSpeed3] = randomTextures2Rotations(opt);
opt.omega   = pi/300;
[ImgSeq4 maxSpeed4] = randomTextureRotation(opt);
opt.expRate = 1.005;
[ImgSeq5 maxSpeed5] = randomTextureExpansion(opt);
rmpath('./../');

% % Display the image sequence.
% for iFrame = 1:size(ImgSeq3,3),
%     imshow(squeeze(ImgSeq3(:,:,iFrame))); drawnow;
%     pause(0.1);
% end

% Display the estimated optic flow.
h       = size(ImgSeq1,1);
w       = size(ImgSeq1,2);
[Y X]   = ndgrid(1:h, 1:w); % pixel coordinates.
sample  = 8;
IndexX  = 1:sample:w;
IndexY  = 1:sample:h; 
% For the display the flow is scaled by division with its maximum speed and
% multiplication with the sampling factor.
len1    = sample/maxSpeed1;
len2    = sample/maxSpeed2;
len3    = sample/maxSpeed3;
len4    = sample/maxSpeed4;
len5    = sample/maxSpeed5;

% *************************************************************************
% Estimate and display flow for the 1st sequence.
% *************************************************************************
[Dx1 Dy1] = estimateOpticFlow2D(ImgSeq1);
figure('Position',[50 50 600 600]);
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx1(IndexY,IndexX,1)*len1, Dy1(IndexY,IndexX,1)*len1,0,'-k');
hold on;
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx1(IndexY,IndexX,2)*len1, Dy1(IndexY,IndexX,2)*len1,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Linear motions.\n',...
    'Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));
drawnow;


% *************************************************************************
% Estimate and display flow for the 2nd sequence.
% *************************************************************************
[Dx2 Dy2] = estimateOpticFlow2D(ImgSeq2);
figure('Position',[150 50 600 600]);
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx2(IndexY,IndexX,1)*len2, Dy2(IndexY,IndexX,1)*len2,0,'-k');
hold on;
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx2(IndexY,IndexX,2)*len2, Dy2(IndexY,IndexX,2)*len2,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Expansion and contraction.\n',...
    'Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));
drawnow;


% *************************************************************************
% Estimate and display flow for the 3rd sequence.
% *************************************************************************
[Dx3 Dy3] = estimateOpticFlow2D(ImgSeq3);
figure('Position',[350 50 600 600]);
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx3(IndexY,IndexX,1)*len3, Dy3(IndexY,IndexX,1)*len3,0,'-k');
hold on;
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx3(IndexY,IndexX,2)*len3, Dy3(IndexY,IndexX,2)*len3,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Clockwise and counterclockwise rotation.\n',...
    'Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));
drawnow;

% *************************************************************************
% Estimate and display flow for the 4th sequence.
% *************************************************************************
[Dx4 Dy4] = estimateOpticFlow2D(ImgSeq4);
figure('Position',[450 50 600 600]);
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx4(IndexY,IndexX,1)*len4, Dy4(IndexY,IndexX,1)*len4,0,'-k');
hold on;
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx4(IndexY,IndexX,2)*len4, Dy4(IndexY,IndexX,2)*len4,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Counterclockwise rotation.\n',...
    'Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));
drawnow;

% *************************************************************************
% Estimate and display flow for the 5th sequence.
% *************************************************************************
[Dx5 Dy5] = estimateOpticFlow2D(ImgSeq5);
figure('Position',[550 50 600 600]);
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx5(IndexY,IndexX,1)*len5, Dy5(IndexY,IndexX,1)*len5,0,'-k');
hold on;
quiver(X(IndexY,IndexX),          Y(IndexY,IndexX),...
       Dx5(IndexY,IndexX,2)*len5, Dy5(IndexY,IndexX,2)*len5,0,'-k');
axis equal ij; axis([-10 w+10 -10 h+10]);
title(sprintf(['Expansion.\n',...
    'Esimated optic flow, sampled %d times for a %d x %d',...
    'pixels resolution.'], sample,h,w));

