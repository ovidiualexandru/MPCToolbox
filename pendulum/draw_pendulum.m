function [J, shapeInserter] = draw_pendulum(height, width, x0, y0, ...
    theta, br, l, t_height, t_length)
%DRAW_PENDULUM Draw a simple pendulum representation.
%   [J, shapeInserter] = draw_pendulum(height, width, x0, y0, theta, ...
%       br, l, t_height, t_length) draws the pendulum and returns the
%       greyscale image in J and the shapeInserter Object in shapeInserter.
%       The image has height in <height>, width in <width>, the cart
%       (trolley) position in (x0,y0), pendulum angle in theta, the ball
%       radius in br, the bar length in l and the cart (trolley) dimensions
%       in (t_height, t_length).
%   Arguments:
%   - height, width: image dimensions
%   - x0, y0: the trolley( cart) position - the base of the pendulum
%   - theta: the pendulum angle relative to the vertical axis.
%   - br: the radius of the ball at the end of the pendulum.
%   - l: the pendulum length
%   - t_height, t_length: the card (trolley) dimensions.
%   Output arguments:
%   - J: the greyscale image with the pendulum
%   - shapeInserter: the vision.ShapeInserter instance used.
%
%   Notes: Computer vision toolbox must be installed.

%% Draw the ball
y1 = y0 - round(cos(theta)*l);
x1 = x0 - round(sin(theta)*l);
circles = [x1 y1 br];
shapeInserter = vision.ShapeInserter('Shape', 'Circles', 'Fill', 1, ...
    'Opacity', 1, 'Antialiasing', 1);
I = ones(height, width); % create a blank image
J = step(shapeInserter, I, circles);
release(shapeInserter);
%% Draw the beam
set(shapeInserter, 'Shape', 'Lines');
J = step(shapeInserter, J, [x0 y0 x1 y1]);
release(shapeInserter);
%% Draw trolley
xtl = x0 - t_length/2;
% xth = x0 + t_length/2;
ytl = y0;
% yth = y0 - t_height;
set(shapeInserter, 'Shape', 'Rectangles', 'Fill', 0);
J = step(shapeInserter, J, [xtl ytl t_length t_height]);
