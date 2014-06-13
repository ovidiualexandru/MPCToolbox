function [J, shapeInserter] = draw_pendulum(height, width, x0, y0, theta, br, l, t_height, t_length)
%% Draw the ball
y1 = y0 - round(cos(theta)*l);
x1 = x0 - round(sin(theta)*l);
circles = [x1 y1 br];
shapeInserter = vision.ShapeInserter('Shape', 'Circles', 'Fill', 1, ...
    'Opacity', 1);
I = ones(height, width); % create a blank image
J = step(shapeInserter, I, circles);
release(shapeInserter);
%% Draw the beam
% shapeInserter = vision.ShapeInserter('Shape', 'Lines', 'Antialiasing', 1);
set(shapeInserter, 'Shape', 'Lines');
J = step(shapeInserter, J, [x0 y0 x1 y1]);
release(shapeInserter);
%% Draw trolley
xtl = x0 - t_length/2;
% xth = x0 + t_length/2;
ytl = y0;
% yth = y0 - t_height;
% shapeInserter = vision.ShapeInserter('Shape', 'Rectangles', 'Antialiasing', 1);
set(shapeInserter, 'Shape', 'Rectangles', 'Fill', 0);
J = step(shapeInserter, J, [xtl ytl t_length t_height]);
