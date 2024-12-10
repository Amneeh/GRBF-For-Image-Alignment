% Visualization and Image Saving for SSD Cost Surface with Varying Sigmas
% This section visualizes the SSD cost surface for each Gaussian sigma value, 
% highlights the estimated minima and true translation, and saves the topview plots 
% as images in the specified output folder.


% --- Visualization and Image Saving ---
output_folder = 'output_images1'; % Folder to save images
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Display results for each sigma
for sigma_idx = 1:length(sigma_values)
    sigma = sigma_values(sigma_idx);
    z_values = z_values_all(:, :, sigma_idx);

    % Find minima for the current sigma (Corrected)
    [~, min_idx] = min(z_values(:));
    [min_i, min_j] = ind2sub(size(z_values), min_idx);

    % Create figure for the surface plot
    figure('Visible', 'off'); % Create figure without displaying
    surf(theta_x_range, theta_y_range, z_values, ...
         'EdgeColor', 'interp', ... % Add interpolated grid lines
         'FaceAlpha', 0.8, ...      % Slight transparency for better visibility
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

    % Set top view
    view(2); % Top-down view
    colormap%('hsv'); % Set colormap for better visualization
    colorbar;

    % Add labels and title
    xlabel('\theta_x (Translation X)');
    ylabel('\theta_y (Translation Y)');
    title(sprintf('Top View SSD with Contours (sigma = %.2f)', sigma));

    % Set figure size and save
    set(gcf, 'Position', [100, 100, 800, 800]); % Set uniform size
    saveas(gcf, fullfile(output_folder, sprintf('delta_0.01_sigma_%.2f_topview_with_lines.png', sigma)));

    % Close figure to save memory
    close(gcf);
end
