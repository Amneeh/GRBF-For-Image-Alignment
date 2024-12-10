% TRANSLATE FUNCTION
% This function applies translation, scaling, rotation, and a transformation matrix 
% to an input image. It returns the transformed image with optional resizing and 
% interpolation for smooth pixel mapping.

function Translated_image = Translate(Og_image, b, scale, theta, c, Z, W, H)
    arguments
        Og_image
        b = [0, 0]
        scale = 1
        theta = 0
        c = [1 0; 0 1]
        Z = 1
        W = []
        H = []
    end

    % Normalize translation vector
    b1 = b(1) / Z;
    b2 = b(2) / Z;

    % Ensure input image is grayscale
    if size(Og_image, 3) == 3
        Og_image = rgb2gray(Og_image);
    end

    % Resize the image
    if ~isempty(W) && ~isempty(H)
        Og_image = imresize(Og_image, [H, W]);
    else
        [H, W] = size(Og_image);
    end

    % Transformation matrix
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)]; % Rotation
    S = scale * eye(2); % Scaling
    T =c; % Combined transformation

    % Center of image
    x0 = (W + 1) / 2;
    y0 = (H + 1) / 2;

    % Initialize output image
    Translated_image = zeros(H, W);

    % Perform the transformation
    for xi2 = 1:H
        for xi1 = 1:W
            % Map pixel coordinates back to original space
            coords = T * ([xi1; xi2] - [x0; y0] - [b1; b2]) + [x0; y0];
            xr = coords(1);
            yr = coords(2);

            % Interpolate pixel values
            if xr >= 1 && xr <= W && yr >= 1 && yr <= H
                Translated_image(xi2, xi1) = interp2(Og_image, xr, yr, 'linear', 0);
            end
        end
    end
end
