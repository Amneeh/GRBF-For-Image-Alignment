% TRANSLATE FUNCTION
% This function applies a combination of transformations, including translation, scaling,
% rotation, and perspective distortion, to an input image. It returns the transformed image 
% and the normalized, resized version of the original image. Optional parameters allow for 
% resizing and depth adjustment during transformation.

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