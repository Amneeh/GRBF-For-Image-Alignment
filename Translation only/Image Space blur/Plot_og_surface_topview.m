% 3D SSD Cost Surface Visualization with Translation Estimation
% This script calculates the Sum of Squared Differences (SSD) cost surface 
% for image registration by evaluating translations of a current image 
% relative to a desired image. It identifies the estimated and true translation 
% points, visualizes the top view of the cost surface, and saves the plot as an image.

% Clean and clear workspace
clc; clear;

%% Define constants
delta = 1; % GRBF spread
n = 2;
colormap('jet'); % Replace 'jet' with your desired colormap

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

plot_original_surface(desiredImage, currentImage, theta_x_range, theta_y_range, output_folder);

function plot_original_surface(desiredImage, currentImage, theta_x_range, theta_y_range, output_folder)
    % Compute the original SSD-based cost function surface
    z_values = zeros(length(theta_y_range), length(theta_x_range)); % Initialize surface

    % Loop through all translation values
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

    % Plot the surface with lines
    figure('Visible', 'off'); % Create figure without displaying
    surf(theta_x_range, theta_y_range, z_values, ...
         'EdgeColor', 'interp', ... % Add interpolated grid lines
         'FaceAlpha', 0.8);         % Slight transparency for better visibility
    hold on;

    % Add contour lines at the base
    contour3(theta_x_range, theta_y_range, z_values, 20, 'k', 'LineWidth', 1);

    % Mark the minima with a red dot
    plot3(minima_x, minima_y, minValue, 'ro', 'MarkerSize', 8, 'LineWidth', 2);

    % Add labels and title
    xlabel('\theta_x (Translation in X)');
    ylabel('\theta_y (Translation in Y)');
    zlabel('Cost Function z(\theta)');
    title('Original Cost Function Surface (SSD)');
    colorbar;
    grid on;
    hold off;

    % Set top view
    view(2); % Top-down view

    % Set figure size
    set(gcf, 'Position', [100, 100, 800, 800]); % Set uniform size

    % Capture the figure as an image
    frame = getframe(gcf);
    img = frame2im(frame);

    % Save the image as JPEG
    save_path = fullfile(output_folder, 'Original_Cost_Function_Surface.jpg');
    imwrite(img, save_path, 'jpg', 'Quality', 95);

    % Close the figure
    close(gcf);
end
