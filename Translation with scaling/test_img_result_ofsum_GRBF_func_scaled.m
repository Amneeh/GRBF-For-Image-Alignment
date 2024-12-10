% IMAGE TRANSFORMATION AND RECOVERY USING GRBF
% This script generates synthetic images by applying translation and scaling transformations 
% to a desired image. It then computes a recovered image using Gaussian Radial Basis Functions (GRBF)
% and the correct values of the transformation. 
% The original, transformed, and recovered images are visualized for comparison.


clc;
clear;
close all
%% Parameters
delta = 0.1; % GRBF spread
sigma = 0.01; % Smoothing kernel variance
n = 1; % Dimension
scale_factor = 0.1; % Scale factor to resize images

translationX = 70; % True translation in X
translationY = 0 * scale_factor; % True translation in Y
translationVector = [translationX, translationY]; % Translation vector
b = -translationVector'; % Inverse translation for testing

%% Generate Synthetic Images
% Load the desired image and preprocess
desiredImage = imread('img1.ppm');
desiredImage = imresize(desiredImage, scale_factor);
true_scale = 2;
a = [1/true_scale 0 ;0 1/true_scale]; % Scaling factors for GRBF
% Apply translation to generate the current image
[desiredImage, currentImage] = Translate(desiredImage, translationVector,true_scale);

%%
Recovered_Image = sum_GRBF(currentImage, delta,b,sigma,n,a); % GRBF representation of desired image

%% Visualize Results

% Visualize the desired and current images
figure;
subplot(1, 3, 1);
imshow(desiredImage, []);
title('Original Image (Normalized)');

subplot(1, 3, 2);
imshow(currentImage, []);
title('Transformed Image');

subplot(1, 3, 3);
imshow(Recovered_Image, []);
title('Recovered Image');



