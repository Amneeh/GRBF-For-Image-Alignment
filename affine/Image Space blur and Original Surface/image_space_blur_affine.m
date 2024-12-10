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

%% Sigma Values to Iterate Over
sigma_values = [0.005, 0.1, 5,15,35]; % Gaussian kernel values

% Create output folder for images
output_folder = 'output_shear_sigma'; % Folder to save images
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% Iterate Over Sigma Values
for sigma = sigma_values
    % Apply Gaussian filtering
    I_c_sigma = imgaussfilt(currentImage, sigma); % Gaussian blur on current image
    I_d_sigma = imgaussfilt(I_d, sigma); % Gaussian blur on desired image

    %% Compute SSD Surface
    z_values = zeros(size(shear_x_grid)); % Initialize SSD surface

    % Initialize waitbar
    total_iterations = numel(shear_x_grid); % Total number of iterations
    current_iteration = 0; % Track the current iteration
    h = waitbar(0, sprintf('Calculating SSD surface for sigma = %.4f...', sigma)); % Initialize waitbar

    for i = 1:numel(shear_x_grid)
        % Current shear values
        shear_x = shear_x_grid(i);
        shear_y = shear_y_grid(i);

        % Create shear transformation matrix
        shear_matrix = [1, shear_x; shear_y, 1];
        tform = affine2d([shear_matrix [0; 0]; 0 0 1]);

        % Apply shear transformation to the desired image
        transformedImage = imwarp(I_d_sigma, tform, 'OutputView', imref2d(size(I_d_sigma)));

        % Compute SSD
        diff = (double(I_c_sigma) - double(transformedImage)).^2;
        z_values(i) = 0.5 * sum(diff(:)); % Sum of Squared Differences

        % Update waitbar
        current_iteration = current_iteration + 1;
        waitbar(current_iteration / total_iterations, h, ...
            sprintf('Calculating SSD surface for sigma = %.4f... %.2f%%', sigma, (current_iteration / total_iterations) * 100));
    end

    % Close waitbar
    close(h);

    %% Find Minima
    [minSSD, minIndex] = min(z_values(:)); % Minimum SSD and its index
    [minRow, minCol] = ind2sub(size(z_values), minIndex); % Row and column indices of minima
    estimated_shear_x = shear_x_range(minCol); % Estimated shear in x-direction
    estimated_shear_y = shear_y_range(minRow); % Estimated shear in y-direction

    fprintf('Sigma = %.4f: True Shear X = %.2f, True Shear Y = %.2f\n', sigma, true_shear_x, true_shear_y);
    fprintf('Sigma = %.4f: Estimated Shear X = %.2f, Estimated Shear Y = %.2f\n', sigma, estimated_shear_x, estimated_shear_y);

    %% Visualization and Image Saving
    % Create 3D surface plot
    figure('Visible', 'off'); % Create figure without displaying
    surf(shear_x_range, shear_y_range, z_values, ...
         'EdgeColor', 'interp', 'FaceAlpha', 0.8, 'FaceColor', 'interp');
    hold on;

    % Compute minima z-values
    true_min_z = z_values(find(shear_y_range == true_shear_y, 1), find(shear_x_range == true_shear_x, 1));
    est_min_z = z_values(minRow, minCol);

    % Add contour lines at the base
    contour3(shear_x_range, shear_y_range, z_values, 20, 'k', 'LineWidth', 1);

    % Mark true minima (green dot)
    plot3(true_shear_x, true_shear_y, true_min_z + 0.1 * abs(true_min_z), ...
          'go', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'g');

    % Mark estimated minima (red dot)
    plot3(estimated_shear_x, estimated_shear_y, est_min_z + 0.1 * abs(est_min_z), ...
          'ro', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    hold off;

    % Set 3D perspective
    view(-45, 30); % Set azimuth and elevation for 3D perspective

    % Add labels and title
    xlabel('Shear X');
    ylabel('Shear Y');
    zlabel('SSD Cost Function');
    title(sprintf('3D SSD Cost Surface for Shear (Sigma = %.4f)', sigma));

    % Add color bar
    colorbar;

    % Set figure size and save
    set(gcf, 'Position', [100, 100, 800, 800]); % Consistent figure size
    saveas(gcf, fullfile(output_folder, sprintf('SSD_3D_Surface_Shear_Sigma_%.4f.png', sigma)));
    close(gcf);

    % Create top-view plot
    figure('Visible', 'off');
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
    title(sprintf('Top View SSD Surface for Shear (Sigma = %.4f)', sigma));
    colorbar;

    % Set figure size and save
    set(gcf, 'Position', [100, 100, 800, 800]); % Consistent figure size
    saveas(gcf, fullfile(output_folder, sprintf('SSD_Top_View_Shear_Sigma_%.4f.png', sigma)));
    close(gcf);
end
