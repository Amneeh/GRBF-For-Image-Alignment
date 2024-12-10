% SHEAR TRANSFORMATION AND SSD VISUALIZATION
% This script computes and visualizes the SSD cost surface for image
% alignment under varying shear transformations with no smoothing. 
% It identifies and compares the true and estimated shear values and saves
% 3D and top-view plots along with visual comparisons of transformed images.

% Clean and clear workspace
clc; clear; close all;

%% Load and preprocess images
I_d = imread('images/nighthawks.jpg'); % Desired image
I_d = rgb2gray(I_d); % Convert to grayscale
scale_factor = 0.1; % Scale factor to resize
I_d = imresize(I_d, scale_factor); % Resize for faster computation

% Parameters
true_shear_x = 0.25; % True shear in x-direction
true_shear_y = -0.3; % True shear in y-direction

% Apply true shear transformation to create a synthetic current image
true_shear_matrix = [1, true_shear_x; true_shear_y, 1]; % Shear transformation matrix
tform = affine2d([true_shear_matrix [0; 0]; 0 0 1]); % Affine2D transformation
currentImage = imwarp(I_d, tform, 'OutputView', imref2d(size(I_d))); % Apply shear

%% Define Shear Ranges
shear_x_range = -0.5:0.01:0.5; % Shear range for x-direction
shear_y_range = -0.5:0.01:0.5; % Shear range for y-direction

% Grid for shear values
[shear_x_grid, shear_y_grid] = meshgrid(shear_x_range, shear_y_range);

%% Compute SSD Surface
z_values = zeros(size(shear_x_grid)); % Initialize SSD surface

% Initialize waitbar
total_iterations = numel(shear_x_grid); % Total number of iterations
current_iteration = 0; % Track the current iteration
h = waitbar(0, 'Calculating SSD surface...'); % Initialize waitbar

for i = 1:numel(shear_x_grid)
    % Current shear values
    shear_x = shear_x_grid(i);
    shear_y = shear_y_grid(i);

    % Create shear transformation matrix
    shear_matrix = [1, shear_x; shear_y, 1];
    tform = affine2d([shear_matrix [0; 0]; 0 0 1]);

    % Apply shear transformation to the desired image
    transformedImage = imwarp(I_d, tform, 'OutputView', imref2d(size(I_d)));

    % Compute SSD
    diff = (double(currentImage) - double(transformedImage)).^2;
    z_values(i) = 0.5 * sum(diff(:)); % Sum of Squared Differences

    % Update waitbar
    current_iteration = current_iteration + 1;
    waitbar(current_iteration / total_iterations, h, ...
        sprintf('Calculating SSD surface... %.2f%%', (current_iteration / total_iterations) * 100));
end

% Close waitbar
close(h);

%% Find Minima
[minSSD, minIndex] = min(z_values(:)); % Minimum SSD and its index
[minRow, minCol] = ind2sub(size(z_values), minIndex); % Row and column indices of minima
estimated_shear_x = shear_x_range(minCol); % Estimated shear in x-direction
estimated_shear_y = shear_y_range(minRow); % Estimated shear in y-direction

fprintf('True Shear X: %.2f, True Shear Y: %.2f\n', true_shear_x, true_shear_y);
fprintf('Estimated Shear X: %.2f, Estimated Shear Y: %.2f\n', estimated_shear_x, estimated_shear_y);

%% --- Visualization and Image Saving ---
output_folder = 'output_shear_images'; % Folder to save images
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Define figure size
figure_size = [100, 100, 800, 800]; % Consistent figure size

% Create 3D surface plot
figure('Visible', 'off'); % Create figure without displaying
surf(shear_x_range, shear_y_range, z_values, ...
     'EdgeColor', 'interp', 'FaceAlpha', 0.8, 'FaceColor', 'interp'); % Smooth surface
hold on;

% Add contour lines at the base
contour3(shear_x_range, shear_y_range, z_values, 20, 'k', 'LineWidth', 1);

% Raise the minima markers slightly above the surface for visibility
true_min_z = z_values(find(shear_y_range == true_shear_y, 1), find(shear_x_range == true_shear_x, 1));
est_min_z = z_values(minRow, minCol);

% Mark true minima (green dot)
plot3(true_shear_x, true_shear_y, true_min_z + 0.1 * abs(true_min_z), ...
      'go', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'g', 'DisplayName', 'True Minima');

% Mark estimated minima (red dot)
plot3(estimated_shear_x, estimated_shear_y, est_min_z + 0.1 * abs(est_min_z), ...
      'ro', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'Estimated Minima');
hold off;

% Set 3D perspective
view(-45, 30); % Set azimuth and elevation for 3D perspective

% Add labels and title
xlabel('Shear X');
ylabel('Shear Y');
zlabel('SSD Cost Function');
title('3D SSD Cost Surface for Shear Transformations');

% Add color bar
colorbar;

% Set figure size and save
set(gcf, 'Position', figure_size); % Consistent figure size
saveas(gcf, fullfile(output_folder, 'SSD_3D_Surface_Shear.png')); % Save as PNG
close(gcf);

% Create top-view plot
figure('Visible', 'off'); % Create figure without displaying
contourf(shear_x_range, shear_y_range, z_values, 20); % Contour plot
hold on;

% Mark true minima (green dot)
plot(true_shear_x, true_shear_y, ...
     'go', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'g');

% Mark estimated minima (red dot)
plot(estimated_shear_x, estimated_shear_y, ...
     'ro', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold off;

% Add labels and title
xlabel('Shear X');
ylabel('Shear Y');
title('Top View SSD Surface for Shear Transformations');
colorbar;

% Set figure size and save
set(gcf, 'Position', figure_size); % Consistent figure size
saveas(gcf, fullfile(output_folder, 'SSD_Top_View_Shear.png')); % Save as PNG
close(gcf);


%%
%% Visualize Transformed Images
% Transform the desired image with the true shear
true_shear_matrix = [1, true_shear_x; true_shear_y, 1];
tform_true = affine2d([true_shear_matrix [0; 0]; 0 0 1]);
transformedImage_true = imwarp(I_d, tform_true, 'OutputView', imref2d(size(I_d)));

% Transform the desired image with the estimated shear
estimated_shear_matrix = [1, estimated_shear_x; estimated_shear_y, 1];
tform_estimated = affine2d([estimated_shear_matrix [0; 0]; 0 0 1]);
transformedImage_estimated = imwarp(I_d, tform_estimated, 'OutputView', imref2d(size(I_d)));

% Visualize the images
figure;
subplot(1, 3, 1);
imshow(I_d, []);
title('Original Image');

subplot(1, 3, 2);
imshow(transformedImage_true, []);
title('True Shear Transformation');

subplot(1, 3, 3);
imshow(transformedImage_estimated, []);
title('Estimated Shear Transformation');

% Save the visualizations
saveas(gcf, fullfile(output_folder, 'Transformed_Images_Comparison.png'));
close(gcf);
