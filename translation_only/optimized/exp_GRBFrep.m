function [newImage, newImage_squared] = exp_GRBFrep(I, delta)
    [height, width] = size(I);
    [X, Y] = meshgrid(1:width, 1:height);
    % Initialize resulting images
    newImage = zeros(height, width);
    newImage_squared = zeros(height, width);
    for i = 1:height
        for j = 1:width
            amplitude = double(I(i, j));
            G = exp(-((X - j).^2 + (Y - i).^2) / (2 * delta^2));
            % Add contributions for I
            newImage = newImage + amplitude * G;
            % Add contributions for I^2
            newImage_squared = newImage_squared + (amplitude^2) * G;
        end
    end
end