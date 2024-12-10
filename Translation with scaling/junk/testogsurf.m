% SSD SURFACE COMPUTATION AND VISUALIZATION
% This script computes and visualizes the Sum of Squared Differences (SSD) cost surface
% for image registration involving translation and scaling. The current image is 
% generated from a desired image using known transformations. The SSD surface is computed 
% over a range of translation and scaling parameters, and the estimated minima are marked 
% on a 3D plot for comparison with the true transformation values.

% Clean and clear workspace
clc; clear; close all;

%% Load and preprocess images
I_d = imread('nighthawks.jpg'); % Desired image
I_d = rgb2gray(I_d); % Convert to grayscale
scale_factor = 0.1; % Scale factor to resize
I_d = imresize(I_d, scale_factor); % Resize for faster computation

% Parameters
translationX = -30 * scale_factor;  % True translation in X
true_scale = 2; % True scale

% Apply true transformation to create a synthetic current image
true_translation_vector = [translationX, 0]; % Ground truth translation
currentImage = imresize(I_d, true_scale); % Apply scaling
currentImage = imtranslate(currentImage, true_translation_vector, 'FillValues', 0); % Apply translation

% Pad current image to match desired image size
[H, W] = size(I_d); % Dimensions of desired image
currentImage = padarray(currentImage, [max(0, H - size(currentImage, 1)), max(0, W - size(currentImage, 2))], 'post');
currentImage = currentImage(1:H, 1:W); % Crop to desired dimensions

%% Define Translation and Scaling Ranges
theta_x_range = -(W - 1)/2:0.1:(W - 1)/2; % Finer translation steps
scale_range = 0.1:0.01:2.5; % Finer scale steps

%% Compute SSD Surface
z_values = zeros(length(scale_range), length(theta_x_range)); % Initialize SSD surface
for i = 1:length(theta_x_range)
    for j = 1:length(scale_range)
        % Current transformation parameters
        theta_x = theta_x_range(i);
        scale = scale_range(j);

        % Apply scaling and translation
        transformedImage = imresize(currentImage, scale); % Apply scaling
        transformedImage = imtranslate(transformedImage, [theta_x, 0], 'FillValues', 0); % Apply translation

        % Pad and crop to match desired image size
        transformedImage = padarray(transformedImage, [max(0, H - size(transformedImage, 1)), max(0, W - size(transformedImage, 2))], 'post');
        transformedImage = transformedImage(1:H, 1:W); % Crop to desired dimensions

        % Compute SSD
        diff = (double(transformedImage) - double(I_d)).^2; % Element-wise difference
        z_values(j, i) = 0.5 * sum(diff(:)); % Sum of Squared Differences
    end
end

%% Find Minima
[minSSD, minIndex] = min(z_values(:)); % Minimum SSD and its index
[minRow, minCol] = ind2sub(size(z_values), minIndex); % Row and column indices of minima
estimated_translation = theta_x_range(minCol); % Estimated translation
estimated_scale = scale_range(minRow); % Estimated scale

fprintf('True Translation: %.2f, True Scale: %.2f\n', translationX, true_scale);
fprintf('Estimated Translation: %.2f, Estimated Scale: %.2f\n', estimated_translation, estimated_scale);

%% Visualize SSD Surface
% 3D Surface Plot
figure;
surf(theta_x_range, scale_range, z_values, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
xlabel('Translation (\theta_x)');
ylabel('Scale (a)');
zlabel('SSD Cost Function');
title('SSD Surface for Translation and Scaling');
colorbar;
view(3);

% Mark the estimated minima
hold on;
plot3(estimated_translation, estimated_scale, minSSD, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
legend('SSD Surface', 'Estimated Minima');
