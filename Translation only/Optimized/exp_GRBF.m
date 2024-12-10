

function [mapping_I_c, mapping_I_c_squared] = exp_GRBF(I_c, delta, b, s, n)
    arguments
        I_c                     % Input image
        delta                   % GRBF spread
        b                       % Translation vector
        s                       % Gaussian spread
        n                       % Dimension parameter
    end

    [height, width] = size(I_c);
    % Scaling factor for the Gaussian
    scaling_factor = (delta / sqrt(delta^2 + s^2))^n;
    % Meshgrid of pixel coordinates
    [X, Y] = meshgrid(1:width, 1:height);
    % Initialize the mappings
    mapping_I_c = zeros(height, width);
    mapping_I_c_squared = zeros(height, width);
    % Ensure image is a double
    if ~isa(I_c, 'double')
        I_c = double(I_c);
    end
    % Vectorize computation for each pixel in the output
    for i = 1:height
        for j = 1:width
            % Compute the translated position
            x_k = [X(i, j); Y(i, j)];
            tau = x_k + b; % Apply scaling and translation
            % Compute distances from tau to all pixels
            distances = (X - tau(1)).^2 + (Y - tau(2)).^2;
            % Gaussian weights
            weights = exp(-distances / (2 * (delta^2 + s^2)));
            % Compute contributions for I_c
            contributions_I_c = I_c .* weights;
            % Compute contributions for I_c.^2
            contributions_I_c_squared = (I_c.^2) .* weights;
            % Accumulate contributions with scaling factor
            mapping_I_c(i, j) = scaling_factor * sum(contributions_I_c(:));
            mapping_I_c_squared(i, j) = scaling_factor * sum(contributions_I_c_squared(:));
        end
    end
end
