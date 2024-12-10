% sum_GRBF FUNCTION
% This function computes a mapping using the Sum of Gaussian Radial Basis Functions (GRBFs). 
% Each pixel in the input image contributes to the mapping (of tau) through a Gaussian 
% kernel centered at a translated position. 
% The output is a smooth representation incorporating translation, scaling, 
% and Gaussian spread. 


function mapping = sum_GRBF(I_c, delta, b, s, n, a)
    arguments
        I_c                     % Input image
        delta                   % GRBF spread
        b                       % Translation vector
        s                       % Gaussian base spread
        n                       % Dimension parameter
        a = [1 0;0 1]              % Optional scaling factors for x and y
    end

    [height, width] = size(I_c);
    % Meshgrid of pixel coordinates
    [X, Y] = meshgrid(1:width, 1:height);
    % Initialize the mapping
    mapping = zeros(height, width);
    % Center the pixel coordinates
    center_x = width / 2;  % Center of the width
    center_y = height / 2; % Center of the height

    % Ensure the input image is a double
    if ~isa(I_c, 'double')
        I_c = double(I_c);
    end

    % Vectorize computation for each pixel in the output
    for i = 1:height
        for j = 1:width
            % Compute the translated position
            x_k = [X(i, j); Y(i, j)];

             tau = a * (x_k - [center_x; center_y]) + [center_x; center_y] + a * b;
       

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
