% COMPUTATION AND VISUALIZATION OF SSD COST FUNCTION WITH GAUSSIAN BLUR
% This script calculates the Sum of Squared Differences (SSD) cost function
% for image alignment under different Gaussian smoothing (sigma) applied to the image space.
% It iterates over translation and scaling parameters, computes the SSD surface,
% and visualizes the results as 3D plots and contour maps, highlighting the 
% true and estimated minima for each sigma value.

% Clean and clear workspace
clc; clear; close all;

%% Load and preprocess images
I_d = imread('nighthawks.jpg'); % Desired image
I_d = rgb2gray(I_d); % Convert to grayscale
scale_factor = 0.1; % Scale factor to resize
I_d = imresize(I_d, scale_factor); % Resize for faster computation

% Parameters
translationX = 30 * scale_factor;  % True translation in X
true_scale = 0.5; % True scale

% Apply true transformation to create a synthetic current image
true_translation_vector = [translationX, 0]; % Ground truth translation
currentImage = imresize(I_d, true_scale); % Apply scaling
currentImage = imtranslate(currentImage, true_translation_vector, 'FillValues', 0); % Apply translation

% Pad current image to match desired image size
[H, W] = size(I_d); % Dimensions of desired image
currentImage = padarray(currentImage, [max(0, H - size(currentImage, 1)), max(0, W - size(currentImage, 2))], 'post');
currentImage = currentImage(1:H, 1:W); % Crop to desired dimensions

%% Define Translation and Scaling Ranges
theta_x_range = -(W - 1)/2:0.5:(W - 1)/2; % Finer translation steps
scale_range = 0.1:0.01:2.5; % Finer scale steps

%% Sigma Values to Iterate Over
sigma_values = [0.1 1 5]; % Gaussian kernel values

% Create output folder for images
output_folder = 'output_images_sigma'; % Folder to save images
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% Iterate Over Sigma Values
for sigma = sigma_values
    % Apply Gaussian filtering
    I_c_sigma = imgaussfilt(currentImage, sigma); % Gaussian blur on current image
    I_d_sigma = imgaussfilt(I_d, sigma); % Gaussian blur on desired image

    %% Compute SSD Surface
    z_values = zeros(length(scale_range), length(theta_x_range)); % Initialize SSD surface

    % Initialize waitbar
    total_iterations = length(scale_range) * length(theta_x_range); % Total number of iterations
    current_iteration = 0; % Track the current iteration
    h = waitbar(0, sprintf('Calculating SSD surface for sigma = %.4f...', sigma)); % Initialize waitbar

    for i = 1:length(theta_x_range)
        for j = 1:length(scale_range)
            % Current transformation parameters
            theta_x = theta_x_range(i);
            scale = scale_range(j);

            % Apply scaling and translation
            transformedImage = imresize(I_d_sigma, scale); % Scale desiredImage
            transformedImage = imtranslate(transformedImage, [theta_x, 0], 'FillValues', 0); % Translate desiredImage

            % Pad and crop to match currentImage size
            transformedImage = padarray(transformedImage, [max(0, H - size(transformedImage, 1)), max(0, W - size(transformedImage, 2))], 'post');
            transformedImage = transformedImage(1:H, 1:W); % Crop to match dimensions

            % Compute SSD
            diff = (double(I_c_sigma) - double(transformedImage)).^2; % Compare against currentImage
            z_values(j, i) = 0.5 * sum(diff(:)); % Sum of Squared Differences

            % Update waitbar
            current_iteration = current_iteration + 1;
            waitbar(current_iteration / total_iterations, h, ...
                sprintf('Calculating SSD surface for sigma = %.4f... %.2f%%', sigma, (current_iteration / total_iterations) * 100));
        end
    end

    % Close waitbar
    close(h);

    %% Find Minima
    [minSSD, minIndex] = min(z_values(:)); % Minimum SSD and its index
    [minRow, minCol] = ind2sub(size(z_values), minIndex); % Row and column indices of minima
    estimated_translation = theta_x_range(minCol); % Estimated translation
    estimated_scale = scale_range(minRow); % Estimated scale

    fprintf('Sigma = %.4f: True Translation = %.2f, True Scale = %.2f\n', sigma, translationX, true_scale);
    fprintf('Sigma = %.4f: Estimated Translation = %.2f, Estimated Scale = %.2f\n', sigma, estimated_translation, estimated_scale);

    %% Visualization and Image Saving
    % Create 3D surface plot (Opposite Perspective)
    figure('Visible', 'off'); % Create figure without displaying
    surf(theta_x_range, scale_range, z_values, ...
         'EdgeColor', 'interp', 'FaceAlpha', 0.8, 'FaceColor', 'interp');
    hold on;
    true_min_z = z_values(find(scale_range == true_scale, 1), find(theta_x_range == translationX, 1));
    est_min_z = z_values(minRow, minCol);

    % Add contour lines at the base
    contour3(theta_x_range, scale_range, z_values, 20, 'k', 'LineWidth', 1);

    % Mark true minima (green dot)
    plot3(translationX, true_scale, true_min_z + 0.1 * abs(true_min_z), ...
          'go', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'g');

    % Mark estimated minima (red dot)
    plot3(estimated_translation, estimated_scale, est_min_z + 0.1 * abs(est_min_z), ...
          'ro', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    hold off;

    % Set 3D perspective (Opposite Side)
    view(-45, 30); % Opposite perspective with azimuth = 135° and elevation = 30°

    % Add labels and title
    xlabel('\theta_x (Translation X)');
    ylabel('Scale (a)');
    zlabel('SSD Cost Function');
    title(sprintf('3D SSD Cost Surface (Opposite Perspective) for Sigma = %.4f', sigma));

    % Add color bar
    colorbar;

    % Set figure size and save
    set(gcf, 'Position', [100, 100, 800, 800]); % Consistent figure size
    saveas(gcf, fullfile(output_folder, sprintf('SSD_3D_Surface_Sigma_%.4f_Opposite.png', sigma))); % Save as PNG
    close(gcf);

    % Create top-view plot
    figure('Visible', 'off');
    contourf(theta_x_range, scale_range, z_values, 20); % Contour plot
    hold on;

    % Mark true minima (green dot)
    plot(translationX, true_scale, ...
         'go', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'g');

    % Mark estimated minima (red dot)
    plot(estimated_translation, estimated_scale, ...
         'ro', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    hold off;

    % Add labels and title
    xlabel('\theta_x (Translation X)');
    ylabel('Scale (a)');
    title(sprintf('Top View SSD Surface for Sigma = %.4f', sigma));
    colorbar;

    % Set figure size and save
    set(gcf, 'Position', [100, 100, 800, 800]); % Consistent figure size
    saveas(gcf, fullfile(output_folder, sprintf('SSD_Top_View_Sigma_%.4f.png', sigma))); % Save as PNG
    close(gcf);
end
