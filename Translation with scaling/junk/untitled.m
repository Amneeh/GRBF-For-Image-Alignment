% Clean and clear workspace
clc; clear;

%% Define constants
delta = 0.005; % GRBF spread
sigma = 0.05; % Fixed sigma value
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

% Initialize cost grid
z_values = zeros(size(theta_x_grid)); % Store cost values for all combinations

%% Compute cost for all (x-translation, scale) combinations
h = waitbar(0, 'Calculating cost grid...');
total_steps = numel(theta_x_grid);

for i = 1:total_steps
    % Current x-translation and scale
    theta_x = theta_x_grid(i); % X-translation
    current_scale = scale_grid(i); % Scale

    % Compute the translation vector and scaling factors
    b = [theta_x; 0]; % Only x-translation, y-translation = 0
    a = [current_scale 0; 0 current_scale]; % Uniform scaling

    % Compute terms of the cost function
    term2 = sum_GRBF(I_c, delta, b, sigma, n, a).* I_d_GRBF;
    term1 = sum_GRBF(I_c.^2, delta, b, sigma, n, a);
    % Compute the total cost for the current combination
    z_values(i) = 0.5 * sum(term1(:)) - sum(term2(:)) + (sqrt(2 * pi * sigma^2) / (2 * (2 * pi * sigma^2)^(n/2)))  * I_d_squared_GRBF_sum;
                  
    % Update waitbar
    if mod(i, round(total_steps / 100)) == 0 % Update every 1%
        waitbar(i / total_steps, h);
    end
end
close(h);

%% Reshape cost grid for visualization
z_values = reshape(z_values, size(theta_x_grid));

%% Visualization
% Plot the 2D surface of the cost function
figure;
surf(theta_x_range, scale_range, z_values, ...
     'EdgeColor', 'none', 'FaceAlpha', 0.8);
xlabel('X Translation (\theta_x)');
ylabel('Scale (a)');
zlabel('Cost Function');
title('Cost Function for X Translation and Scale');
colorbar;
view(3);

% Plot the contour map of the cost function
figure;
contourf(theta_x_range, scale_range, z_values, 20);
xlabel('X Translation (\theta_x)');
ylabel('Scale (a)');
title('Contour Plot of Cost Function');
colorbar;

% Highlight true translation and scale
hold on;
plot(-translationX, 1/true_scale, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'True Minima');
legend('Contours', 'True Minima');
hold off;

%% Find and visualize the minimum
% Locate the minimum cost
[min_cost, min_idx] = min(z_values(:));
[min_row, min_col] = ind2sub(size(z_values), min_idx);

% Extract the parameters corresponding to the minimum cost
estimated_theta_x = theta_x_range(min_col); % X-translation
estimated_scale = scale_range(min_row); % Scale

% Apply the estimated parameters to the desired image
[~, estimated_image] = Translate(I_c, [estimated_theta_x, 0], estimated_scale);

%% Visualization of the estimated minima
figure;

% Display the current image
subplot(1, 3, 1);
imshow(I_c, []);
title('Current Image');

% Display the desired image
subplot(1, 3, 2);
imshow(I_d, []);
title('Desired Image');

% Display the estimated transformed image
subplot(1, 3, 3);
imshow(estimated_image, []);
title(sprintf('Transformed Image (\\theta_x = %.2f, Scale = %.2f)', ...
    estimated_theta_x, estimated_scale));

% Print the estimated parameters and cost
fprintf('Minimum Cost: %.4f\n', min_cost);
fprintf('Estimated Translation (X): %.2f\n', estimated_theta_x);
fprintf('Estimated Scale: %.2f\n', estimated_scale);
