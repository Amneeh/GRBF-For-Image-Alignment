% sum_GRBF FUNCTION
% This function computes a mapping using the Sum of Gaussian Radial Basis Functions (GRBFs) representation. 
% Each pixel in the input image contributes to the mapping (tau) through a Gaussian 
% kernel centered at a translated position, scaled by a specified factor. 
% The output is a smooth representation incorporating translation, scaling, 
% and Gaussian spread.
function mapping = sum_GRBF(I_c, delta, b, s, n, a, c)
    arguments
        I_c                     % Input image
        delta                   % GRBF spread
        b                       % Translation vector
        s                       % Gaussian base spread
        n                       % Dimension parameter
        a = [1; 1]              % Optional scaling factors for x and y
        c = [1 0; 0 1]          % Shear matrix
    end

    [height, width] = size(I_c);
    [X, Y] = meshgrid(1:width, 1:height);
    mapping = zeros(height, width);

    if ~isa(I_c, 'double')
        I_c = double(I_c);
    end
    center_x = width / 2;  % Center of the width
    center_y = height / 2; % Center of the height

    a = diag(a);
    for i = 1:height
        for j = 1:width
            x_k = [X(i, j); Y(i, j)];

            %%%%%%%%%%%%% Problem %%%%%%%%%%%%%%%
            tau = c * (x_k - [center_x; center_y]) + [center_x; center_y] + c .* b;

            % Compute distances as 2D matrix
            distances = (X - tau(1)).^2 + (Y - tau(2)).^2;

            % Compute sigma^2
            sigma_squared = s^2 * (1 + norm(x_k)^2);

            % Update scaling factor
            scaling_factor = (delta / sqrt(delta^2 + sigma_squared))^n;

            % Compute weights as 2D matrix
            weights = exp(-distances / (2 * sigma_squared));

            % Compute contributions as 2D matrix
            contributions = I_c .* weights;

            % Accumulate contributions
            mapping(i, j) = sum(contributions(:) .* scaling_factor(:));
        end
    end
end
