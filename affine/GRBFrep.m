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
