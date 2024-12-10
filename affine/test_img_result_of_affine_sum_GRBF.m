clc;
clear;
%% 

% Mean Squared Error: 1.5009
% PSNR Value: -1.7635 dB
% Normalized Cross-Correlation (NCC) Value: 0.18863
% >> SHx 0.6
% Mean Squared Error: 0.82145
% PSNR Value: 0.85419 dB
% Normalized Cross-Correlation (NCC) Value: 0.53434
% >> SHx 0.25
% Mean Squared Error: 1.2844
% PSNR Value: -1.087 dB
% Normalized Cross-Correlation (NCC) Value: 0.25815
% >> SHx 0
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


%%
function Translated_image = Translate(Og_image, b, scale, theta, c, Z, W, H)
    arguments
        Og_image
        b = [0, 0]
        scale = 1
        theta = 0
        c = [1 0; 0 1]
        Z = 1
        W = []
        H = []
    end

    % Normalize translation vector
    b1 = b(1) / Z;
    b2 = b(2) / Z;

    % Ensure input image is grayscale
    if size(Og_image, 3) == 3
        Og_image = rgb2gray(Og_image);
    end

    % Resize the image
    if ~isempty(W) && ~isempty(H)
        Og_image = imresize(Og_image, [H, W]);
    else
        [H, W] = size(Og_image);
    end

    % Transformation matrix
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)]; % Rotation
    S = scale * eye(2); % Scaling
    T = S * R * c; % Combined transformation

    % Center of image
    x0 = (W + 1) / 2;
    y0 = (H + 1) / 2;

    % Initialize output image
    Translated_image = zeros(H, W);

    % Perform the transformation
    for xi2 = 1:H
        for xi1 = 1:W
            % Map pixel coordinates back to original space
            coords = T * ([xi1; xi2] - [x0; y0] - [b1; b2]) + [x0; y0];
            xr = coords(1);
            yr = coords(2);

            % Interpolate pixel values
            if xr >= 1 && xr <= W && yr >= 1 && yr <= H
                Translated_image(xi2, xi1) = interp2(Og_image, xr, yr, 'linear', 0);
            end
        end
    end
end
