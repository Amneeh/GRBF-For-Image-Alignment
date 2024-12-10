% COMPUTATION AND VISUALIZATION OF COST FUNCTION USING GRBF FOR TRANSLATION AND SCALING
% This script calculates the cost function for image transformation under 
% different Gaussian Radial Basis Function (GRBF) smoothing values (sigma). 
% It evaluates the cost function for combinations of x-translation and scale 
% transformations and visualizes the results as 3D surfaces and contour plots. 
% True and estimated transformation values are marked for comparison.


% Clean and clear workspace
clc; clear;

%% Define constants
delta = 0.001; % GRBF spread
sigma_values = [0.001 0.005 0.01 0.05]; % Different sigma values to iterate over
n = 2; % Dimensions

% Load and preprocess images
I_d = imread('nighthawks.jpg'); % Desired image
I_d = rgb2gray(I_d);
scale_factor = 0.05; % Scale factor to resize
I_d = imresize(I_d, scale_factor);
translationX = -30 * scale_factor;  % True translation in X
true_scale = 2; % True scale 

% Apply known translation and scale
[I_d, I_c] = Translate(I_d, [translationX, 0], true_scale); % Apply translation in X only
[H, W] = size(I_d); % Image dimensions

% Define x-translation and scale ranges
theta_x_range = -(W - 1)/2:2:(W - 1)/2; % Translation range in X
scale_range = 0.1:0.05:2.5; % Scale range
[theta_x_grid, scale_grid] = meshgrid(theta_x_range, scale_range);
% Precompute GRBF Representations of the desired image
I_d_GRBF = GRBFrep(I_d, delta); % GRBF representation of desired image
I_d_squared_GRBF = GRBFrep(I_d.^2, delta); % GRBF representation of squared desired image
I_d_squared_GRBF_sum = sum(I_d_squared_GRBF(:)); % Precomputed scalar

% Initialize cost grids
z_values_all = []; % Store cost values for all sigma values

%% Iterate over sigma values
for sigma_idx = 1:length(sigma_values)
    sigma = sigma_values(sigma_idx); % Current sigma
    disp(['Processing sigma = ', num2str(sigma)]);

    % Initialize cost grid for current sigma
    z_values = zeros(size(theta_x_grid));

    %% Compute cost for all (x-translation, scale) combinations
    h = waitbar(0, ['Calculating cost grid for sigma = ', num2str(sigma)]);
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

    % Reshape and store cost grid for current sigma
    z_values = reshape(z_values, size(theta_x_grid));
    z_values_all = cat(3, z_values_all, z_values); % Store all results
end

%% Visualization for each sigma
for sigma_idx = 1:length(sigma_values)
    sigma = sigma_values(sigma_idx);
    z_values = z_values_all(:, :, sigma_idx);

    % Plot the 2D surface
    figure;
    surf(theta_x_range, scale_range, z_values, ...
         'EdgeColor', 'none', 'FaceAlpha', 0.8);
    xlabel('X Translation (\theta_x)');
    ylabel('Scale (a)');
    zlabel('Cost Function');
    title(['Cost Function for Sigma = ', num2str(sigma)]);
    colorbar;
    view(3);

    % Plot the contour map
    figure;
    contourf(theta_x_range, scale_range, z_values, 20);
    xlabel('X Translation (\theta_x)');
    ylabel('Scale (a)');
    title(['Contour Plot of Cost Function for Sigma = ', num2str(sigma)]);
    colorbar;

    % Locate the minimum for current sigma
    [min_cost, min_idx] = min(z_values(:));
    [min_row, min_col] = ind2sub(size(z_values), min_idx);
    estimated_theta_x = theta_x_range(min_col); % X-translation
    estimated_scale = scale_range(min_row); % Scale

    % Highlight true and estimated values on contour plot
    hold on;
    plot(-translationX, 1/true_scale, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'True Minima');
    plot(theta_x_range(min_col), scale_range(min_row), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Estimated Minima');
    legend('Contours', 'True Minima', 'Estimated Minima');
    hold off;

    % Display results for this sigma
    fprintf('Sigma = %.2f: Min Cost = %.4f, Estimated X = %.2f, Estimated Scale = %.2f\n', ...
        sigma, min_cost, estimated_theta_x, estimated_scale);
end
