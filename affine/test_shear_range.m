clc;
clear;
close all;

%% Parameters
delta = 0.01;                    % GRBF spread
sigma_values = [0.01]; % Array of sigma values to evaluate smoothing
n = 2;                        % Dimension
scale_factor = 0.1;           % Scale factor to resize images

% True translation (set to zero for shear testing)
translationX = 0;
translationY = 0;
translationVector = [translationX, translationY]; % Translation vector
b = -translationVector'; % Inverse translation for testing

% Shear ranges
shear_values = 0:0.1:1; % Shear values for c12 and c21
[shear_x_grid, shear_y_grid] = meshgrid(shear_values, shear_values);

%% Generate Synthetic Images
% Load the desired image and preprocess
desiredImage = imread('images/nighthawks.jpg');
desiredImage = imresize(desiredImage, scale_factor);

true_scale = 1;               % True scale
a = [1/true_scale; 1/true_scale]; % Scaling factors for GRBF

% Apply a known shear to generate the current image
true_shear = [1, 0.25; -0.3, 1]; % Known true shear matrix
[desiredImage, currentImage] = Translate(desiredImage, translationVector, true_scale, 0, true_shear);

figure;
% Desired and Current Images
subplot(2, 1, 1);
imshow(desiredImage, []);
title('Desired Image');

subplot(2, 1, 2);
imshow(currentImage, []);
title('Current Image');
