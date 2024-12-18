% SSD SURFACE COMPUTATION AND VISUALIZATION
% This script computes the Sum of Squared Differences (SSD) cost surface
% for image alignemnt over a range of shearing in the x and y direction.
% It generates a transformed image by applying known translation and scaling
% to a desired image, computes the SSD surface, visualizes it using 3D and
% contour plots, and highlights the true and estimated minima.
clc;
clear;
close all;

%% Parameters
delta = 1;                    % GRBF spread
sigma_values = [0.005 0.009 0.01 0.02 ]; % Array of sigma values to evaluate smoothing
n = 2;                        % Dimension
scale_factor = 0.1;           % Scale factor to resize images

% True translation (set to zero for shear testing)
translationX = 0;
translationY = 0;
translationVector = [translationX, translationY]; % Translation vector
b = -translationVector'; % Inverse translation for testing

% Shear ranges
shear_values = -0.9:0.01:0.9; % Shear values for c12 and c21
[shear_x_grid, shear_y_grid] = meshgrid(shear_values, shear_values);

%% Generate Synthetic Images
% Load the desired image and preprocess
desiredImage = imresize(double(rgb2gray(imread('images/nighthawks.jpg'))), scale_factor) / 255;
true_scale = 1;               % True scale
a = [1/true_scale; 1/true_scale]; % Scaling factors for GRBF

% Apply a known shear to generate the current image
true_shear = [1, 0.25; -0.3, 1]; % Known true shear matrix
currentImage = Translate(desiredImage, translationVector, true_scale, 0, true_shear);

%% Precompute GRBF Representations of the desired image
desiredImage_GRBF = GRBFrep(desiredImage, delta);
I_d_squared_GRBF = GRBFrep(desiredImage.^2, delta); % GRBF representation of squared desired image
I_d_squared_GRBF_sum = sum(I_d_squared_GRBF(:)); % Precomputed scalar


%% Loop Through Sigma Values
cost_grid_all = []; % Store cost grids for all sigma values
for sigma_idx = 1:length(sigma_values)
    sigma = sigma_values(sigma_idx);
    cost_grid = zeros(size(shear_x_grid)); % Initialize cost grid
        h = waitbar(0, sprintf('Processing sigma = %.2f', sigma));

    % Loop through all shear combinations
    for i = 1:numel(shear_x_grid)
        c12 = shear_x_grid(i); % Shear in X
        c21 = shear_y_grid(i); % Shear in Y
        shear_matrix = [1, c12; c21, 1];
        shear_inverse = inv(shear_matrix);

        % Compute the cost value for the current shear
        term1 = sum_GRBF(currentImage.^2, delta, b, sigma, n, a, shear_inverse);
        term2 = sum_GRBF(currentImage, delta, b, sigma, n, a, shear_inverse) .* desiredImage_GRBF;

        cost_grid(i) = 0.5 * sum(term1(:)) - sum(term2(:))+ (sqrt(2 * pi * sigma^2) / (2 * (2 * pi * sigma^2)^(n/2)))  * I_d_squared_GRBF_sum;;
                waitbar(i / numel(shear_x_grid), h, ...
            sprintf('Processing sigma = %.2f (%d/%d)', sigma, i, numel(shear_x_grid)));

    end

    % Store the computed cost grid
    cost_grid = reshape(cost_grid, size(shear_x_grid));
    cost_grid_all = cat(3, cost_grid_all, cost_grid);
end

%% Visualization
true_shear_x = 0.25; % True shear in X direction (c12)
true_shear_y = -0.3; % True shear in Y direction (c21)

for sigma_idx = 1:length(sigma_values)
    sigma = sigma_values(sigma_idx);
    cost_grid = cost_grid_all(:, :, sigma_idx);

    % Find the estimated minima (location of the minimum cost value)
    [min_val, min_idx] = min(cost_grid(:)); % Find the minimum value and index
    [min_row, min_col] = ind2sub(size(cost_grid), min_idx); % Convert index to row and column
    estimated_shear_x = shear_x_grid(1, min_col); % Estimated shear in X
    estimated_shear_y = shear_y_grid(min_row, 1); % Estimated shear in Y

    % Find the closest indices for the true shear values
    [~, true_shear_x_idx] = min(abs(shear_values - true_shear_x)); % Closest index for true c12
    [~, true_shear_y_idx] = min(abs(shear_values - true_shear_y)); % Closest index for true c21

    % 3D Surface Plot
    figure;
    surf(shear_x_grid, shear_y_grid, cost_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    hold on;

    % Plot the estimated minima (red marker)
    plot3(estimated_shear_x, estimated_shear_y, min_val, ...
          'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Estimated Minima');

    % Plot the true shear values (green marker)
    true_val = cost_grid(true_shear_y_idx, true_shear_x_idx); % True value at the true shear indices
    plot3(true_shear_x, true_shear_y, true_val, ...
          'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'True Shear');

    hold off;

    % Labels and Title
    xlabel('Shear in X Direction (c_{12})');
    ylabel('Shear in Y Direction (c_{21})');
    zlabel('Cost Function');
    title(sprintf('Cost Function (sigma = %.2f)', sigma));
    colorbar;
    legend('Cost Surface', 'Estimated Minima', 'True Shear');
    view(3);

    % Contour Plot
    figure;
    contourf(shear_x_grid, shear_y_grid, cost_grid, 20);
    hold on;

    % Plot the estimated minima (red marker)
    plot(estimated_shear_x, estimated_shear_y, ...
         'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Estimated Minima');

    % Plot the true shear values (green marker)
    plot(true_shear_x, true_shear_y, ...
         'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'True Shear');

    hold off;

    % Labels and Title
    xlabel('Shear in X Direction (c_{12})');
    ylabel('Shear in Y Direction (c_{21})');
    title(sprintf('Cost Function Contour (sigma = %.2f)', sigma));
    colorbar;
    legend('Contours', 'Estimated Minima', 'True Shear');
end

% Load the desired image
desiredImage = imread('images/nighthawks.jpg'); % Adjust path as needed
desiredImage = rgb2gray(imresize(desiredImage, scale_factor));

% Apply the true shear transformation
true_shear = [1, 0.25; -0.3, 1]; % True shear matrix
true_sheared_image = sum_GRBF(currentImage,delta,-[0, 0]',sigma,2,[1;1],inv(true_shear));

[min_val, min_idx] = min(cost_grid(:));
[min_row, min_col] = ind2sub(size(cost_grid), min_idx);
estimated_shear_x = shear_x_grid(1, min_col);
estimated_shear_y = shear_y_grid(min_row, 1);

% Apply the estimated shear transformation (identity matrix for now)
estimated_shear = [1, estimated_shear_x; estimated_shear_y, 1]; % Replace with the estimated shear matrix
estimated_sheared_image= sum_GRBF(currentImage,delta,-[0, 0]',sigma,2,[1;1],inv(estimated_shear));

% Visualize the images
figure;

subplot(2, 2, 1);
imshow(desiredImage, []);
title('Original Desired Image');

subplot(2, 2, 2);
imshow(currentImage, []);
title('Current Image');

subplot(2, 2, 3);
imshow(true_sheared_image, []);
title('True Sheared Image');

subplot(2, 2, 4);
imshow(estimated_sheared_image, []);
title('Estimated Sheared Image');
