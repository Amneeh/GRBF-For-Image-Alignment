clc;
clear;
profile on 
%% Parameters
delta = 1; % GRBF spread
sigma = 0.001; % Smoothing kernel variance
n = 2; % Dimension
scale_factor = 0.6; % Scale factor to resize images

translationX = 0 * scale_factor; % True translation in X
translationY = 0 * scale_factor; % True translation in Y
translationVector = [translationX, translationY]; % Translation vector
b = -translationVector'; % Inverse translation for testing

%% Generate Synthetic Images
% Load the desired image and preprocess
desiredImage = imread('images/apple.jpg');
desiredImage = imresize(desiredImage, scale_factor);
true_scale = 1;
a = [1/true_scale;1/true_scale]; % Scaling factors for GRBF
% Apply translation to generate the current image
% Apply shear and its inverse
shear = [1, 0.1; 0.1, 1];
[desiredImage1, currentImage1] = Translate(desiredImage, translationVector, true_scale, 0, shear);
[desiredImage, currentImage] = Translate(currentImage1, translationVector, true_scale, 0, inv(shear));

figure;

% Visualize the desired and current images
subplot(1, 4, 1);
imshow(desiredImage1, []);
title('Desired Image');

subplot(1, 4, 2);
imshow(currentImage1, []);
title('Current Image');

subplot(1, 4, 3);
imshow(desiredImage, []);
title('Desired Image');

subplot(1, 4, 4);
imshow(currentImage, []);
title('Current Image');

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
        c = [1 0; 0 1]              % Shear matrix
        Z = 1                       % Depth 
        W = []                      % Optional new width for resizing
        H = []                      % Optional new height for resizing
    end

    % Preprocess translation parameters
    b1 = b(1) / Z;
    b2 = b(2) / Z;

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

    % Perform translation, scaling, rotation, and shear
    for xi1 = 1:W
        for xi2 = 1:H
            % Compute the transformed coordinates
            transformed_coords = c * ([(xi1 - x0); (xi2 - y0)]) + [x0; y0] + [b1; b2];
            xr = transformed_coords(1);
            yr = transformed_coords(2);

            % Interpolate pixel value at (xr, yr)
            if xr >= 1 && xr <= W && yr >= 1 && yr <= H
                Translated_image(xi2, xi1) = interp2(Desired_image, xr, yr, 'bilinear'); 
            end
        end
    end
end

profile viewer