% Visualization and Image Saving for SSD Cost Surface with Varying Sigmas
% This section visualizes the SSD cost surface for each Gaussian sigma value, 
% highlights the estimated minima and true translation, and saves the 3D plots 
% as images in the specified output folder.

% --- Visualization and Image Saving ---
output_folder = 'output_images1'; % Folder to save images
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Define consistent figure size
figure_size = [100, 100, 800, 800]; % Consistent figure size (width x height)

% Display results for each sigma
for sigma_idx = 1:length(sigma_values)
    sigma = sigma_values(sigma_idx);
    z_values = z_values_all(:, :, sigma_idx);

    % Find minima for the current sigma
    [~, min_idx] = min(z_values(:));
    [min_i, min_j] = ind2sub(size(z_values), min_idx);

    % Create figure for the surface plot
    figure('Visible', 'off'); % Create figure without displaying
    surf(theta_x_range, theta_y_range, z_values, ...
         'EdgeColor', 'interp', ... % Add interpolated grid lines
         'FaceAlpha', 0.9, ...      % Slight transparency for better visibility
         'FaceColor', 'interp');    % Smooth face colors
    hold on;

    % Add contour lines to the surface
    contour3(theta_x_range, theta_y_range, z_values, 20, 'k', 'LineWidth', 1);

    % Mark the estimated minima (red dot)
    plot3(theta_x_range(min_j), theta_y_range(min_i), z_values(min_i, min_j), ...
          'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

    % Mark the true translation point (green dot)
    plot3(theta_x_range(true_i), theta_y_range(true_j), z_values(true_j, true_i), ...
          'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    hold off;

    % Set flipped 3D perspective
    view(-45, 30); % Negative azimuth to flip perspective (adjust as needed)

    % Add labels, title
    xlabel('\theta_x (Translation X)');
    ylabel('\theta_y (Translation Y)');
    zlabel('SSD Cost Function');
    title(sprintf('3D SSD Cost Surface (sigma = %.2f)', sigma));

    % Remove color bar
    % colorbar; <-- Removed for a clean plot

    % Set figure size and save
    set(gcf, 'Position', figure_size); % Set uniform figure size
    saveas(gcf, fullfile(output_folder, sprintf('delta_0.01_sigma_%.2f_3d.png', sigma)));

    % Close figure to save memory
    close(gcf);
end
