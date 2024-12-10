% Visualization of GRBF-Based Image represenation of a current image using
% the Sum_GRBF function.
% This script generates an image pair using translation, calculates
% GRBF-based representations, and visualizes the results. It highlights 
% differences between the desired and current image.
% To confirm the functions work properly, the resutling image from Term2 with the
% correct translation should match the desired image with an addition of a blurring affect.

clc;
clear;

%% Parameters
delta = 1; % GRBF spread
sigma = 0.5; % Smoothing kernel variance
n = 2; % Dimension

scale_factor = 0.1; % Scale factor to resize

translationX = 150 * scale_factor;  % True translation in X
translationY = 200 * scale_factor; % True translation in Y
translationVector = [translationX, translationY]; % Translation applied to one image
b = -translationVector'; % Translation for testing

%% Generate Synthetic Images
% Load the desired image and preprocess
desiredImage = imresize(double(rgb2gray(imread('img1.ppm'))), scale_factor) / 255;
% Apply translation to generate the current image
currentImage = Translate(desiredImage,[translationX, translationY]); 
%% Term Calculations

Term2 = sum_GRBF(currentImage, delta, -translationVector', sigma, n); % Second term


%% Visualize Results
figure;
% Desired and Current Images
subplot(2, 2, 1);
imshow(desiredImage, []);
title('Desired Image');

subplot(2, 2, 2);
imshow(currentImage, []);
title('Current Image');

% GRBF Visualizations
subplot(2, 2, 3);
imshow(Term2, []);
title('Eq.20 Representation of Current image');
desiredImage_rescaled = desiredImage * max(Term2(:)) / max(desiredImage(:));
diff = abs(desiredImage_rescaled - Term2);
% GRBF Visualizations
subplot(2, 2, 4);
imshow(diff, []);
title('difference');
