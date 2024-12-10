% IMAGE TRANSFORMATION AND RECOVERY USING GRBF
% This script applies translation and shear transformations to an image 
% and attempts to recover it using Gaussian Radial Basis Functions (GRBF). 
% It visualizes the original, transformed, and recovered images, 
% and computes the difference between the original and recovered images.
clc;
clear;
%% 
delta = 0.75;    
sigma = 0.005;
% Load and preprocess the image
Og_image = imread('images/img1.ppm');
scale_factor = 0.1; % Scale factor to resize
Og_image = imresize(Og_image, scale_factor);
Og_image = double(rgb2gray(Og_image)) / 255;

% Define transformation parameters
scale = 1; % Scaling factor
shear = [1, 0.25 ; -0.3, 1]; % Shear matrix
theta = 0; % No rotation
estimate = [1 0.25; -0.3 1];
c = inv(estimate); % No perspective distortion
b = [100, 10]; % Translation vector

% Apply translation
Transformed_image= Translate(Og_image,b,scale,0,shear);

% Apply reverse translation
% Recovered_image = Translate(Transformed_image,[0 0],scale,0,c);%
Recovered_image = sum_GRBF(Transformed_image,delta,b,sigma,2,[1/scale;1/scale],c);
% Display results
figure;
subplot(1, 3, 1);
imshow(Og_image, []);
title('Original Image (Normalized)');

subplot(1, 3, 2);
imshow(Transformed_image, []);
title('Transformed Image');

subplot(1, 3, 3);
imshow(Recovered_image, []);
title('Recovered Image');
Og_image = Og_image * max(Recovered_image(:)) / max(Og_image(:));
% Compute and visualize difference
diff = abs(Recovered_image - Og_image);
figure;
imshow(diff, []);
title('Difference Image');
% 
% %%
% % Structural Similarity Index
% [ssimval, ssimmap] = ssim(Recovered_image, Og_image);
% figure;
% imshow(ssimmap, []);
% title(['SSIM Map (SSIM Value = ', num2str(ssimval), ')']);
% 
% % Mean Squared Error
% mseval = immse(Recovered_image, Og_image);
% disp(['Mean Squared Error: ', num2str(mseval)]);
% 
% % Peak Signal-to-Noise Ratio
% psnrval = psnr(Recovered_image, Og_image);
% disp(['PSNR Value: ', num2str(psnrval), ' dB']);
% 
% % Normalized Cross-Correlation
% norm_corr = normxcorr2(Og_image, Recovered_image);
% figure;
% imshow(norm_corr, []);
% title('Normalized Cross-Correlation');
% 
% % Visualize Difference as Heatmap
% figure;
% imagesc(diff);
% colorbar;
% title('Difference Heatmap');
% %% NCC 
% % Subtract mean to make images zero mean
% Og_image_zero_mean = Og_image - mean(Og_image(:));
% Recovered_image_zero_mean = Recovered_image - mean(Recovered_image(:));
% 
% % Compute NCC
% ncc_value = sum(Og_image_zero_mean(:) .* Recovered_image_zero_mean(:)) / ...
%             sqrt(sum(Og_image_zero_mean(:).^2) * sum(Recovered_image_zero_mean(:).^2));
% disp(['Normalized Cross-Correlation (NCC) Value: ', num2str(ncc_value)]);
% 
% % Visualize NCC as a map
% ncc_map = normxcorr2(Og_image_zero_mean, Recovered_image_zero_mean);
% figure;
% imagesc(ncc_map);
% colorbar;
% title('Normalized Cross-Correlation Map');
