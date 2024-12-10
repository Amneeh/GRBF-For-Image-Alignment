% SSD SURFACE COMPUTATION AND VISUALIZATION FOR TRANSLATION AND SCALING
% This script computes the Sum of Squared Differences (SSD) cost surface
% for image alignment over a range of x-translation and scale parameters.
% It generates a transformed image by applying known translation and scaling
% to a desired image, computes the SSD surface, visualizes it using 3D and
% contour plots, and highlights the true and estimated minima.


% Clean and clear workspace
clc; clear;

%% Define constants
scale_factor = 0.1; % Scale factor to resize
translationX = -30 * scale_factor; % True translation in X
true_scale = 2; % True scale

% Load and preprocess images
I_d = imread('nighthawks.jpg'); % Desired image
I_d = rgb2gray(I_d); % Convert to grayscale
I_d = imresize(I_d, scale_factor); % Downscale for faster processing
[H, W] = size(I_d); % Image dimensions

% Define x-translation and scale ranges
theta_x_range = -(W - 1)/2:1:(W - 1)/2; % Translation range in X
scale_range = 0.1:0.05:2.5; % Scale range
[theta_x_grid, scale_grid] = meshgrid(theta_x_range, scale_range);

% Apply the known translation and scaling to generate the current image
[I_d, I_c] = Translate(I_d, [translationX, 0], true_scale); % Apply translation in X only
%% Compute SSD surface
z_values = zeros(size(theta_x_grid)); % Initialize cost grid

% Iterate over all translations and scales
total_iterations = length(theta_x_range) * length(scale_range);
current_iteration = 0;
h = waitbar(0, 'Calculating SSD surface...');
for i = 1:numel(theta_x_grid)
    % Current translation and scale
    theta_x = theta_x_grid(i);
    current_scale = scale_grid(i);
    
    % Apply translation and scaling
    [~, I_transformed] = Translate(I_c, [theta_x, 0], current_scale); % Apply known translation and scaling
    % Compute SSD for the current transformation
    diff =( double(I_transformed) - double(I_d)).^2;
    z_values(i) = sum(diff(:));
    current_iteration = current_iteration + 1;
    waitbar(current_iteration / total_iterations, h, sprintf('Calculating SSD surface... %.2f%%', (current_iteration / total_iterations) * 100));
end

% Reshape the cost grid
z_values = reshape(z_values, size(theta_x_grid));

%% Visualize the SSD surface
% 3D Surface Plot
figure;
surf(theta_x_range, scale_range, z_values, ...
     'EdgeColor', 'none', 'FaceAlpha', 0.8); % Smooth surface
xlabel('X Translation (\theta_x)');
ylabel('Scale (a)');
zlabel('SSD Cost Function');
title('SSD Surface');
colorbar;
view(3);

% Contour Plot
figure;
contourf(theta_x_range, scale_range, z_values, 20); % Contour plot
xlabel('X Translation (\theta_x)');
ylabel('Scale (a)');
title('Contour Plot of SSD Surface');
colorbar;

% Locate and highlight the minima
[min_cost, min_idx] = min(z_values(:));
[min_row, min_col] = ind2sub(size(z_values), min_idx);
estimated_theta_x = theta_x_range(min_col);
estimated_scale = scale_range(min_row);

% Highlight true and estimated values
hold on;
plot(-translationX, 1/true_scale, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'True Minima');
plot(estimated_theta_x, estimated_scale, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Estimated Minima');
legend('Contours', 'True Minima', 'Estimated Minima');
hold off;

% Display results
fprintf('True Translation: %.2f, True Scale: %.2f\n', -translationX, 1/true_scale);
fprintf('Estimated Translation: %.2f, Estimated Scale: %.2f\n', estimated_theta_x, estimated_scale);

% Compute the SSD at the estimated minima
[~, I_at_minima] = Translate(I_c, [estimated_theta_x, 0], estimated_scale);

% Plot the images
figure;

% Original Desired Image
subplot(1, 3, 1);
imshow(I_d, []);
title('Desired Image');

% Current Image (Transformed)
subplot(1, 3, 2);
imshow(I_c, []);
title('Current Image (I_c)');

% Image at the Estimated Minima
subplot(1, 3, 3);
imshow(I_at_minima, []);
title(sprintf('Image at Estimated Minima\n(tx = %.2f, scale = %.2f)', ...
    estimated_theta_x, estimated_scale));

% Highlight difference between desired and estimated minima image
figure;
imshowpair(I_d, I_at_minima, 'diff');
title(sprintf('Difference Between Desired and Image at Minima\n(tx = %.2f, scale = %.2f)', ...
    estimated_theta_x, estimated_scale));
colorbar;

%%
disp(class(I_d));          % Check the type of desired image
disp(class(I_at_minima));  % Check the type of transformed image
disp(class(I_c));          % Check the type of the current image