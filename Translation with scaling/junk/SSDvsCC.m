% Clean and clear workspace
clc; clear;

%% Define constants
delta = 0.005; % GRBF spread
sigma = 0.01; % Fixed sigma value
n = 2; % Dimensions

% Load and preprocess images
I_d = imread('nighthawks.jpg'); % Desired image
I_d = rgb2gray(I_d);
scale_factor = 0.1; % Scale factor to resize
I_d = imresize(I_d, scale_factor);
translationX = -30 * scale_factor;  % True translation in X
true_scale = 2; % True scale 

% Apply known translation and scale
[I_d, I_c] = Translate(I_d, [translationX, 0], true_scale); % Apply translation in X only
[H, W] = size(I_d); % Image dimensions

% Define x-translation and scale ranges
theta_x_range = -(W - 1)/4:1:(W - 1)/4; % Translation range in X
scale_range = 0.1:0.05:2.5; % Scale range
[theta_x_grid, scale_grid] = meshgrid(theta_x_range, scale_range);

% Precompute GRBF Representations of the desired image
I_d_GRBF = GRBFrep(I_d, delta); % GRBF representation of desired image
I_d_squared_GRBF = GRBFrep(I_d.^2, delta); % GRBF representation of squared desired image
I_d_squared_GRBF_sum = sum(I_d_squared_GRBF(:)); % Precomputed scalar

% Initialize cost grids
z_values_1 = zeros(size(theta_x_grid)); % Cost grid for cost function 1
z_values_2 = zeros(size(theta_x_grid)); % Cost grid for cost function 2

%% Test Cost Function 1
disp('Testing Cost Function 1...');
tic; % Start timer
for i = 1:numel(theta_x_grid)
    theta_x = theta_x_grid(i); % X-translation
    current_scale = scale_grid(i); % Scale
    b = [theta_x; 0];
    a = [current_scale 0; 0 current_scale];
    
    term2 = sum_GRBF(I_c, delta, b, sigma, n, a);
    z = (term2 - I_d_GRBF).^2;
    z_values_1(i) = sum(z(:));
end
time_cost_function_1 = toc; % End timer
fprintf('Cost Function 1 Time: %.4f seconds\n', time_cost_function_1);

%% Test Cost Function 2
disp('Testing Cost Function 2...');
tic; % Start timer
for i = 1:numel(theta_x_grid)
    theta_x = theta_x_grid(i); % X-translation
    current_scale = scale_grid(i); % Scale
    b = [theta_x; 0];
    a = [current_scale 0; 0 current_scale];
    
    term2 = sum_GRBF(I_c, delta, b, sigma, n, a);
    z = term2 .* I_d_GRBF;
    z_values_2(i) = sum(z(:));
end
time_cost_function_2 = toc; % End timer
fprintf('Cost Function 2 Time: %.4f seconds\n', time_cost_function_2);

%% Reshape cost grids
z_values_1 = reshape(z_values_1, size(theta_x_grid));
z_values_2 = reshape(z_values_2, size(theta_x_grid));

%% Visualization of Cost Function Results
% Plot Cost Function 1
figure;
surf(theta_x_range, scale_range, z_values_1, ...
     'EdgeColor', 'none', 'FaceAlpha', 0.8);
xlabel('X Translation (\theta_x)');
ylabel('Scale (a)');
zlabel('Cost Function 1');
title('Cost Function 1 (SSD-Based)');
colorbar;
view(3);

% Plot Cost Function 2
figure;
surf(theta_x_range, scale_range, -z_values_2, ...
     'EdgeColor', 'none', 'FaceAlpha', 0.8);
xlabel('X Translation (\theta_x)');
ylabel('Scale (a)');
zlabel('Cost Function 2');
title('Cost Function 2 (Product-Based)');
colorbar;
view(3);

%% Locate and Visualize Minima
% Locate minima for both cost functions
[min_cost_1, min_idx_1] = min(z_values_1(:));
[min_row_1, min_col_1] = ind2sub(size(z_values_1), min_idx_1);

[min_cost_2, min_idx_2] = min(-z_values_2(:));
[min_row_2, min_col_2] = ind2sub(size(-z_values_2), min_idx_2);

% Extract the parameters corresponding to the minima
estimated_theta_x_1 = theta_x_range(min_col_1); % X-translation for cost function 1
estimated_scale_1 = scale_range(min_row_1); % Scale for cost function 1

estimated_theta_x_2 = theta_x_range(min_col_2); % X-translation for cost function 2
estimated_scale_2 = scale_range(min_row_2); % Scale for cost function 2

% Display estimated results
fprintf('Cost Function 1 Minima: Translation (X) = %.2f, Scale = %.2f, Cost = %.4f\n', ...
    estimated_theta_x_1, estimated_scale_1, min_cost_1);
fprintf('Cost Function 2 Minima: Translation (X) = %.2f, Scale = %.2f, Cost = %.4f\n', ...
    estimated_theta_x_2, estimated_scale_2, min_cost_2);

% Overlay minima on contour plots
figure;
contourf(theta_x_range, scale_range, z_values_1, 20);
hold on;
plot(estimated_theta_x_1, estimated_scale_1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('X Translation (\theta_x)');
ylabel('Scale (a)');
title('Cost Function 1 Contour with Minima');
colorbar;

figure;
contourf(theta_x_range, scale_range, -z_values_2, 20);
hold on;
plot(estimated_theta_x_2, estimated_scale_2, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
xlabel('X Translation (\theta_x)');
ylabel('Scale (a)');
title('Cost Function 2 Contour with Minima');
colorbar;
