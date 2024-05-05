clear
clc
clf

% Define the overall dimensions of the matrix, assuming m + r = 5 for the given coordinates
m = 3; % m is not explicitly defined, but 5 - r (which is 2) will be used for m
n = 5; % n is defined by the provided coordinates
r = 2; % r is not explicitly defined, but the last coordinate pair suggests it's 2

% Create the figure
figure(1);
set(gcf , 'Color' , 'w')
hold on

% Define the coordinates of the custom polygon vertices
custom_polygon_x = [0, 5, 5, 4, 0, 0];
custom_polygon_y = [5, 5, 0, 0, 4, 5];

% Fill the custom polygon with gray color
fill(custom_polygon_x, custom_polygon_y, [0.8 0.8 0.8], 'EdgeColor', 'none', 'LineWidth', 2);

% four lines
line([0, 0], [0, 5], 'Color', 'k', 'LineWidth', 2);
line([0, 5], [5, 5], 'Color', 'k', 'LineWidth', 2);
line([5, 5], [5, 0], 'Color', 'k', 'LineWidth', 2);
line([5, 0], [0, 0], 'Color', 'k', 'LineWidth', 2);

% Add horizontal line at y = r
line([0, n], [r, r], 'Color', 'k', 'LineWidth', 2);

% Add vertical line at x = r
line([r, r], [0, m+r], 'Color', 'k', 'LineWidth', 2);

% Set the axis limits and aspect ratio
axis equal;
axis off;
xlim([0, n]);
ylim([0, m+r]);

% Font size for the labels
fontSize = 12; % Change this value as needed

% Add text labels for each section E11, E12, E21, and E22
text(r/2, m/2 + r, '$E_{11}$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'latex', 'FontSize', fontSize);
text((n+r)/2, m/2 + r, '$E_{12}$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'latex', 'FontSize', fontSize);
text(r/2, (r)/2, '$E_{21}$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'latex', 'FontSize', fontSize);
text((n+r)/2, (r)/2, '$E_{22}$', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'latex', 'FontSize', fontSize);


% Add a left curly brace
annotation('textbox', [0.17, 0.765, 0.05, 0.165], 'String', '$\left\{\rule{0pt}{0.8cm}\right.$', 'Interpreter', 'latex', 'FontSize', fontSize, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

% Add a right curly brace at the red area
annotation('textbox', [0.84, 0.2, 0.05, 0.6], 'String', '$\left.\rule{0pt}{1.6cm}\right\}$', 'Interpreter', 'latex', 'FontSize', fontSize, 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

% Show the figure
drawnow;
