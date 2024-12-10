% 3D SSD Cost Surface Visualization with Translation Estimation
% This script calculates the Sum of Squared Differences (SSD) cost surface 
% for image registration by evaluating translations of a current image 
% relative to a desired image. It identifies the estimated and true translation 
% points, visualizes the 3D cost surface, and saves the plot as an image.

% Clean and clear workspace
clc; clear;

%% Define constants
delta = 1; % GRBF spread
n = 2;
colormap('jet'); % Use desired colormap

% Load and preprocess images
I_d = imread('nighthawks.jpg'); % Desired image
I_d = rgb2gray(I_d);
scale_factor = 0.1; % Scale factor to resize
I_d = imresize(I_d, scale_factor);
translationX = 75 * scale_factor;  % True translation in X
translationY = -50 * scale_factor; % True translation in Y
translationVector = [translationX, translationY]; % Ground truth translation

[desiredImage, currentImage] = Translate(I_d, translationVector); % Apply known translation

%% Translation Ranges
[H, W] = size(desiredImage); % Image dimensions
theta_x_range = -(W - 1)/2:(W - 1)/2; % Translation range for X
theta_y_range = -(H - 1)/2:(H - 1)/2; % Translation range for Y

%% Plot Original Surface and Save
output_folder = 'output_images'; % Folder to save images
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Compute SSD-based cost function
z_values = zeros(length(theta_y_range), length(theta_x_range)); % Initialize surface

for i = 1:length(theta_x_range)
    for j = 1:length(theta_y_range)
        % Apply translation
        b = [theta_x_range(i), theta_y_range(j)];
        translatedImage = imtranslate(currentImage, -b, 'FillValues', 0); % Align current image

        % Compute SSD cost
        diff = (translatedImage - desiredImage).^2;
        z_values(j, i) = 0.5 * sum(diff(:)); % Sum of Squared Differences
    end
end

% Find the minima of the cost function
[minValue, minIndex] = min(z_values(:)); % Minimum value and index
[minRow, minCol] = ind2sub(size(z_values), minIndex); % Convert index to row and column
minima_x = theta_x_range(minCol); % X-coordinate of the minima
minima_y = theta_y_range(minRow); % Y-coordinate of the minima

% Create a figure for the surface plot
figure('Visible', 'off'); % Create figure without displaying
surf(theta_x_range, theta_y_range, z_values, ...
     'EdgeColor', 'interp', ... % Add interpolated grid lines
     'FaceAlpha', 0.9, ...      % Slight transparency for better visibility
     'FaceColor', 'interp');    % Smooth face colors
hold on;

% Add contour lines to the surface
contour3(theta_x_range, theta_y_range, z_values, 20, 'k', 'LineWidth', 1);

% Mark the estimated minima (red dot)
plot3(minima_x, minima_y, minValue, ...
      'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Mark the true translation point (green dot)
plot3(translationX, translationY, z_values(theta_y_range == translationY, theta_x_range == translationX), ...
      'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');

hold off;

% Add labels and title
xlabel('\theta_x (Translation X)');
ylabel('\theta_y (Translation Y)');
zlabel('SSD Cost Function');
title('3D SSD Cost Surface');

% Set flipped 3D perspective
view(-45, 30); % Negative azimuth to flip perspective (adjust as needed)

% Define consistent figure size
figure_size = [100, 100, 800, 800]; % Consistent figure size (width x height)
set(gcf, 'Position', figure_size); % Set uniform figure size

% Save as an image
save_path = fullfile(output_folder, '3D_Cost_Surface.png');
saveas(gcf, save_path);

% Close the figure to save memory
close(gcf);
