% GRBFrep FUNCTION
% This function computes a Gaussian Radial Basis Function (GRBF) representation
% of an input image. Each pixel in the input image is treated as a center of 
% a Gaussian, weighted by its intensity. The resulting image is the sum of all 
% these Gaussian contributions.

function newImage = GRBFrep(I, delta)
    [height, width] = size(I);
    [X, Y] = meshgrid(1:width, 1:height);
    % initialize resulting image
    newImage = zeros(height, width);
    for i = 1:height
        for j = 1:width
            amplitude = double(I(i, j));
            G = amplitude * exp(-((X - j).^2 + (Y - i).^2) / (2 * delta^2));
            newImage = newImage + G;
        end
    end
end
