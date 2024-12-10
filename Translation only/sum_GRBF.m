% sum_GRBF FUNCTION
% This function computes a mapping using the Sum of Gaussian Radial Basis Functions (GRBFs). 
% Each pixel in the input image contributes to the mapping (of tau) through a Gaussian 
% kernel centered at a translated position. 
% The output is a smooth representation incorporating translation, scaling, 
% and Gaussian spread. 


function mapping = sum_GRBF(I_c, delta, b, s, n)
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
    % Initialize the mapping
    mapping = zeros(height, width);
    %ensure image is a double
    if ~isa(I_c, 'double')
    I_c = double(I_c);
    end
    % Vectorize computation for each pixel in the output
    for i = 1:height
        for j = 1:width
            % % % % % % % % % % % Compute the translated position
            x_k = [X(i, j); Y(i, j)];
            tau = x_k + b; % Apply scaling and translation
            % Compute distances from tau to all pixels
            distances = (X - tau(1)).^2 + (Y - tau(2)).^2;
            % Gaussian weights
            weights = exp(-distances / (2 * (delta^2 + s^2)));
            % Compute contributions
            contributions = I_c .* weights;
            % Accumulate contributions with scaling factor
            mapping(i, j) = scaling_factor * sum(contributions(:));
        end
    end
end
