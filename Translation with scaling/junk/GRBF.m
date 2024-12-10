% Clean and clear workspace
clc; clear;
profile on
%% Define constants
delta = 1; % GRBF spread
n = 2;
% Define the array of sigma values
sigma = [0.01];%,0.5,5,15,35];

%% Load and preprocess images
I_d = imread('Images\apple.jpg'); % Desired image
I_d = rgb2gray(I_d);
scale_factor = 0.3; % Scale factor to resize
I_d = imresize(I_d, scale_factor);
translationX = 18 * scale_factor;  % True translation in X
translationY = -27 * scale_factor; % True translation in Y
true_scale = 2; % True Scale 

%use translate function to apply known translation
[I_d,I_c] = Translate(I_d,[translationX, translationY],true_scale); % Apply known translation
[H, W] = size(I_d); % Image dimensions
theta_x_range = -(W - 1)/2:2:(W - 1)/2; % Translation range (x)
theta_y_range = -(H - 1)/2:2:(H - 1)/2; % Translation range (y)
scale_range = 0.5;%:0.5:1.5;  % Scale rang (a)
[theta_x_grid, theta_y_grid] = meshgrid(theta_x_range, theta_y_range);

% Initialize cost grid
z_values_all = []; % Store cost values for all sigmas

%% Precompute GRBF Representations of the desired image (constant in loop since it is sigma dependant)
I_d_GRBF = GRBFrep(I_d, delta); % GRBF representation of desired image
I_d_squared_GRBF = GRBFrep(I_d.^2, delta); % GRBF representation of squared desired image
I_d_squared_GRBF_sum = sum(I_d_squared_GRBF(:)); % Precomputed scalar


%% Loop through scale values
for scale_idx = 1:length(scale_range)
    % Current scale
    current_scale = scale_range(scale_idx);
    
    % Initialize cost grid for current scale
    z_values = zeros(size(theta_x_grid)); 
    %% Vectorize translation parameters
    b = [theta_x_grid(:)'; theta_y_grid(:)']; % Translation vector (x and y)
    a = [current_scale;current_scale]; % Apply uniform scaling
    % Pre-allocate
    term1 = zeros([size(I_d), size(b, 2)]);
    term2 = zeros([size(I_d), size(b, 2)]);

    % Create waitbar
    h = waitbar(0, sprintf('Calculating cost function for scale = %.2f', current_scale));

    % Loop through all translations
    for k = 1:size(b, 2)
        % Compute terms for the current translation parameter
        term1(:, :, k)  = sum_GRBF(I_c.^2, delta, b(:, k), sigma, n, a);
        term2(:, :, k)  = sum_GRBF(I_c, delta, b(:, k), sigma, n, a) .* I_d_GRBF;
        % Update waitbar
        waitbar(k / size(b, 2), h); 
    end

    close(h); 
    %%
    % Compute the SSD cost for all translations at once
    z_values = 0.5 * sum(term1, [1, 2]) - sum(term2, [1, 2]) + (sqrt(2 * pi * sigma^2) / (2 * (2 * pi * sigma^2)^(n/2)))  * I_d_squared_GRBF_sum;
    z_values = reshape(z_values, size(theta_x_grid));

    % Store z_values for the current sigma
    z_values_all = cat(3, z_values_all, z_values); 
end

%% --- Visualization ---

% Find closest indices to true translation
[~, true_i] = min(abs(theta_x_range + translationX)); % Closest index for X translation
[~, true_j] = min(abs(theta_y_range + translationY)); % Closest index for Y translation

% Ensure valid indices before plotting
if isempty(true_i) || isempty(true_j)
    error('True translation indices could not be matched. Check the ranges.');
end

% Display results for each scale
for scale_idx = 1:length(scale_range)
    current_scale = scale_range(scale_idx);
    z_values = z_values_all(:, :, scale_idx);

    % Find minima for the current scale
    [~, min_idx] = min(z_values(:));
    [min_i, min_j] = ind2sub(size(z_values), min_idx); 

    figure; 
    % Plot the 3D surface
    surf(theta_x_range, theta_y_range, z_values, ...
        'EdgeColor', 'interp', ...
        'FaceAlpha', 0.8);
    hold on;

    % Add contour lines
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
    title(sprintf('3D Surface Plot of SSD (Scale = %.2f)', current_scale));
    view(3);

    % Add a legend for clarity
    legend('Surface with Grid', 'Contours', 'Estimated Translation', 'True Translation', 'Location', 'best');
end

profile viewer


function mapping = sum_GRBF(I_c, delta, b, s, n, a)
    arguments
        I_c                     % Input image
        delta                   % GRBF spread
        b                       % Translation vector
        s                       % Gaussian base spread
        n                       % Dimension parameter
        a = [1; 1]              % Optional scaling factors for x and y
    end

    [height, width] = size(I_c);
    % Meshgrid of pixel coordinates
    [X, Y] = meshgrid(1:width, 1:height);
    % Initialize the mapping
    mapping = zeros(height, width);
    center = [width / 2; height / 2];
    % Ensure the input image is a double
    if ~isa(I_c, 'double')
        I_c = double(I_c);
    end
    % Vectorize computation for each pixel in the output
    for i = 1:height
        for j = 1:width
            % Compute the translated position
            x_k = [X(i, j); Y(i, j)];
            tau = a .* (x_k - center) + center + a .* b;

            % Compute distances from tau to all pixels
            distances = (X - tau(1)).^2 + (Y - tau(2)).^2;

            % Define sigma^2 as a function of each pixel's coordinate
            sigma_squared = s^2 * (1 + (X.^2 + Y.^2));

            % Update scaling factor with variable sigma^2
            scaling_factor = (delta ./ sqrt(delta^2 + sigma_squared)).^n;

            % Gaussian weights with variable sigma^2
            weights = exp(-distances ./ (2 * sigma_squared));

            % Compute contributions
            contributions = I_c .* weights;

            % Accumulate contributions with pixel-specific scaling factor
            mapping(i, j) = sum(contributions(:) .* scaling_factor(:));
        end
    end
end

% FUNCTION TRANSLATE: 
% Takes an image and any optional transformation parameters to return an
% Image translated as per the inputed parameters.
% Function also returns the inputed images normalized and resized. 
% Function can also resize image optionaly if W and H are entered.
function [Desired_image, Translated_image] = Translate(Og_image, b, scale, theta, c, Z, W, H)
    arguments
        Og_image                    % Input image
        b = [0, 0]                  % Translation vector [tx, ty]
        scale = 1                   % Scaling factor 
        theta = 0                   % Rotation angle in radians
        c = [0, 0]                  % Perspective distortion coefficients
        Z = 1                       % Depth 
        W = []                      % Optional new width for resizing
        H = []                      % Optional new height for resizing
    end

    % Preprocess translation parameters with depth adjustment
    b1 = b(1) / Z;
    b2 = b(2) / Z;
    c1 = c(1);
    c2 = c(2);

    % Define transformation matrix A
    a11 = cos(theta) * scale;
    a12 = sin(theta) * scale;
    a21 = -sin(theta) * scale;
    a22 = cos(theta) * scale;
    A = [a11 a12; a21 a22];

    % Ensure the input image is grayscale
    if size(Og_image, 3) == 3
        Og_image = rgb2gray(Og_image);
    end

    % Resize the image if W and H are provided
    if ~isempty(W) && ~isempty(H)
        Og_image = imresize(Og_image, [H, W]);
    else
        [H, W] = size(Og_image);
    end

    % Normalize the original image
    Desired_image = double(Og_image) / 255;

    % Initialize the translated image
    Translated_image = zeros(H, W);

    % Center of the image
    x0 = fix(W / 2);
    y0 = fix(H / 2);

    % Perform translation, scaling, rotation, and perspective distortion
    for xi1 = 1:W
        for xi2 = 1:H
            % Compute the transformed coordinates
            v1 = A(1, 1) * (xi1 - x0) + A(1, 2) * (xi2 - y0) + x0 + b1;
            v2 = A(2, 1) * (xi1 - x0) + A(2, 2) * (xi2 - y0) + y0 + b2;
            v3 = 1 + c1 * (xi1 - x0) + c2 * (xi2 - y0);

            xr = v1 / v3;
            yr = v2 / v3;
    
            % Interpolate pixel value at (xr, yr)
            if xr >= 1 && xr <= W && yr >= 1 && yr <= H
                Translated_image(xi2, xi1) = interp2(Desired_image, xr, yr, 'bilinear'); 
            end
        end
    end
end