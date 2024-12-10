% Image Registration using SSD with Gaussian Filtering in the image space and Visualizations
% This script computes the Sum of Squared Differences (SSD) cost function 
% for image registration under various Gaussian blur scales (sigma) that
% are applied to the image space.
% It identifies the optimal translation parameters that align a 
% reference image with a translated image and visualizes the results.

% Clean and clear workspace
clc; clear;

%% Define constants
delta = 0.0001; % GRBF spread
n = 1; %dimensions
%Define the array of sigma values
sigma_values = [0.01,0.5,5,15,35]; % Example sigma values

% Load and preprocess images
I_d = imread('nighthawks.jpg'); % Desired image
I_d = rgb2gray(I_d);
scale_factor = 0.1; % Scale factor to resize
I_d = imresize(I_d, scale_factor);

translationX = 75 * scale_factor;  % True translation in X
translationY = -50 * scale_factor; % True translation in Y

[I_d,I_c] = Translate(I_d,[translationX, translationY]); % Apply known translation

[H, W] = size(I_d); % Image dimensions
theta_x_range = -(W - 1)/2:(W - 1)/2; % Translation range (x)
theta_y_range = -(H - 1)/2:(H - 1)/2; % Translation range (y)
[theta_x_grid, theta_y_grid] = meshgrid(theta_x_range, theta_y_range);

% Initialize cost grid
z_values_all = []; % Store cost values for all sigmas


%% Precompute images convoluted with a kernel
% Loop through sigma values
for sigma_idx = 1:length(sigma_values)
    sigma = sigma_values(sigma_idx); % Current sigma
    cost_grid = zeros(size(theta_x_grid)); % Initialize cost grid
    I_c_sigma = imgaussfilt(I_c, sigma);
    I_d_sigma = imgaussfilt(I_d, sigma);

    % Create waitbar
    h = waitbar(0, sprintf('Processing sigma = %.2f', sigma));

    for i = 1:length(theta_x_range)
        for j = 1:length(theta_y_range)
            % Apply translation
            b = [theta_x_range(i), theta_y_range(j)];
            I_c_sigma_translated = imtranslate(I_c_sigma, -b, 'FillValues', 0); % Align current image

            % Compute SSD cost
            diff = (I_c_sigma_translated - I_d_sigma).^2;
            cost_grid(j, i) = 0.5 * sum(diff(:)); % Sum of Squared Differences
        end

        % Update waitbar
        waitbar(i / length(theta_x_range), h, ...
            sprintf('Processing sigma = %.2f (%d/%d)', sigma, i, length(theta_x_range)));

    end
    % % Nested loop to iterate over all combinations of theta_x and theta_y
    % for i = 1:length(theta_x_range) % Loop over theta_x
    %     for j = 1:length(theta_y_range) % Loop over theta_y
    %         % Current translation parameters
    %         theta_x = theta_x_range(i); % Translation in X
    %         theta_y = theta_y_range(j); % Translation in Y
    %         translation_vector = [theta_x; theta_y];
    %         I_c_sigma = imtranslate(I_c_sigma,[theta_x, theta_y]);
    % 
    %         % Store cost in the corresponding grid location
    %         ssd = (I_c_sigma - I_d_sigma).^2 * 0.5 ;
    %         cost_grid(j, i) = sum(ssd(:));
    %     end
    % 
    %     % Update waitbar
    %     waitbar(i / length(theta_x_range), h, ...
    %         sprintf('Processing sigma = %.2f (%d/%d)', sigma, i, length(theta_x_range)));
    % end
    % 
    % Close the waitbar
    close(h);

    % Store the computed cost grid for the current sigma
    z_values_all = cat(3, z_values_all, cost_grid);
end

%% --- Visualization ---

% Find closest indices to true translation
[~, true_i] = min(abs(theta_x_range + translationX)); % Closest index for X translation
[~, true_j] = min(abs(theta_y_range + translationY)); % Closest index for Y translation

% Ensure valid indices before plotting
if isempty(true_i) || isempty(true_j)
    error('True translation indices could not be matched. Check the ranges.');
end

% Display results for each sigma
for sigma_idx = 1:length(sigma_values)
    sigma = sigma_values(sigma_idx);
    z_values = z_values_all(:, :, sigma_idx);

    % Find minima for the current sigma (Corrected)
    [~, min_idx] = min(z_values(:));
    [min_i, min_j] = ind2sub(size(z_values), min_idx); 

    figure; 
% Cost Function Visualization with Minima and True Translation
nexttile;

% Plot the 3D surface with enhanced visuals
surf(theta_x_range, theta_y_range, z_values, ...
     'EdgeColor', 'interp', ... % Add interpolated grid lines
     'FaceAlpha', 0.8);         % Slight transparency for better visibility
hold on;

% Add contour lines for better context
contour3(theta_x_range, theta_y_range, z_values, 20, 'k', 'LineWidth', 1);

% Mark the estimated minima (red dot)
plot3(theta_x_range(min_j), theta_y_range(min_i), z_values(min_i, min_j), ...
      'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Estimated Translation');

% Mark the true translation point (green dot)
plot3(theta_x_range(true_i), theta_y_range(true_j), z_values(true_j, true_i), ...
      'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'True Translation');

hold off;

% Add color bar and labels
colorbar;
xlabel('\theta_x (Translation X)');
ylabel('\theta_y (Translation Y)');
zlabel('SSD Cost Function');
title(sprintf('3D Surface Plot of SSD (sigma = %.2f)', sigma));
view(3);

% Add a legend for clarity
legend('Surface with Grid', 'Contours', 'Estimated Translation', 'True Translation', 'Location', 'best');

end
