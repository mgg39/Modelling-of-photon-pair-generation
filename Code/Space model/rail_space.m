% Rail positions and dimensions ----------------------
% currently modelled and rectangles - TODO: dev
rail_dimensions = [10, 1, 1]; % [length, width, height] - 3D space
rail2_dimensions = rail1_dimensions;


rail1_position = [0];
rail2_position = [2];


%Figure ----------------------------------------------

% Create a figure for the 3D plot
figure;
hold on;

% Plot the first rail as a 3D rectangle
drawRect3D(rail1_position, rail1_dimensions, 'b');

% Plot the second rail as a 3D rectangle
drawRect3D(rail2_position, rail2_dimensions, 'r')

% Customize the plot with labels, titles, etc.
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Dual Rail Setup');

% Set axis limits if needed
% axis([-5, 5, -5, 5, 0, 2]);

% Turn off hold to prevent further additions to the plot
hold off;
