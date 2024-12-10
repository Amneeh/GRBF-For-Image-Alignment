% --- Visualization and Image Saving ---
output_folder = 'output_images_shear'; % Folder to save images
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Save true values for visualization
true_shear_x = true_shear_x; % True X translation
true_shear_y = true_shear_y; % True inverse scale

% Display results for each sigma
for sigma_idx = 1:length(sigma_values)
    sigma = sigma_values(sigma_idx);
    z_values = cost_grid_all(:, :, sigma_idx);

    % Find minima for the current sigma
    [~, min_idx] = min(z_values(:)); % Index of the minimum value in the grid
    [min_row, min_col] = ind2sub(size(z_values), min_idx); % Row and column indices of minima
    estimated_shear_x = shear_values(min_col); % X-translation corresponding to the minima
    estimated_shear_y= shear_values(min_row); % Scale corresponding to the minima

    % Debug minima values
    fprintf('Sigma = %.4f: True shear x = %.2f, True shear y = %.2f\n', sigma, true_shear_x, true_shear_y);
    fprintf('Sigma = %.4f: Estimated shear x = %.2f, Estimated shear y = %.2f\n', sigma, estimated_shear_x, estimated_shear_y);

    % Create figure for the surface plot
    figure('Visible', 'off'); % Create figure without displaying
    surf(shear_values, shear_values, z_values, ...
         'EdgeColor', 'interp', ... % Add interpolated grid lines
         'FaceAlpha', 0.8, ...      % Slight transparency for better visibility
         'FaceColor', 'interp');    % Smooth face colors
    hold on;

    % Add contour lines to the surface
    contour3(shear_values, shear_values, z_values, 20, 'k', 'LineWidth', 1);

    % Mark the estimated minima (red dot)
    plot3(estimated_shear_x, estimated_shear_y, z_values(min_row, min_col), ...
          'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Estimated Minima');

    % Mark the true translation point (green dot)
    plot3(true_shear_x, true_shear_y, min(z_values(:)), ...
          'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'True Minima');
    hold off;

    % Set top view
    view(2); % Top-down view
    colormap; % Set colormap for better visualization
    colorbar;

    % Add labels and title
    xlabel('Shear in x ');
    ylabel('Shear in y ');
    title(sprintf('Top View SSD with Contours (sigma = %.4f)', sigma));

    % Set figure size and save
    set(gcf, 'Position', [100, 100, 800, 800]); % Set uniform size
    saveas(gcf, fullfile(output_folder, sprintf('delta_%.4f_sigma_%.4f_topview_with_lines.png', delta, sigma)));

    % Close figure to save memory
    close(gcf);
end

