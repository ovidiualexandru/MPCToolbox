addpath('./pendulum');
addpath('./util');
%% Tunables
vid_width = 1000;
vid_height = 700;
%% Preparations
t = 0:N-1;
figh = figure(2);
set(gcf, 'Units', 'pixels');
set(gcf, 'Position', [100 100 vid_width vid_height]);
%% Draw static images
h1 = subplot(4,2,1);
plot(t,u);
rescaleYLim(gca, [du(2) du(1)]*1.1);
grid on
title('Input u');
line([0;N],[du(1);du(1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[du(2);du(2)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
l1 = line([t(1);t(1)],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

subplot(4,2,2);
plot(t,X(1,:));
rescaleYLim(gca, [dx(2,1) dx(1,1)]*1.1);
grid on
title('Arm position x_1');
line([0;N],[dx(1,1);dx(1,1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[dx(2,1);dx(2,1)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
l2 = line([t(1);t(1)],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

subplot(4,2,4);
plot(t,X(2,:));
rescaleYLim(gca, [dx(2,2) dx(1,2)]*1.1);
grid on
title('Arm speed x_2');
line([0;N],[dx(1,2);dx(1,2)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[dx(2,2);dx(2,2)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
l3 = line([t(1);t(1)],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

subplot(4,2,6);
plot(t,X(3,:));
grid on
rescaleYLim(gca, [dx(2,3) dx(1,3)]*1.1); 
title('Trolley position x_3');
line([0;N],[dx(1,3);dx(1,3)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[dx(2,3);dx(2,3)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
l4 = line([t(1);t(1)],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

subplot(4,2,8);
plot(t,X(4,:));
rescaleYLim(gca, [dx(2,4) dx(1,4)]*1.1);
grid on
title('Trolley speed x_4');
line([0;N],[dx(1,4);dx(1,4)], 'LineStyle', '--', 'Color', [1 0 0]); %%Upper bound
line([0;N],[dx(2,4);dx(2,4)], 'LineStyle', '--', 'Color', [1 0 0]); %%Lower bound
l5 = line([t(1);t(1)],get(gca,'YLim'), 'LineStyle', '--', 'Color', [0 1 0]);

axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
titlestring = sprintf('Run with N=%d, Nc = %d',N,Nc);
text(0.5, 1, ['\bf ' titlestring],'HorizontalAlignment','center','VerticalAlignment', 'top');
%% Draw dynamic images and save to file
subplot(4,2, [5,7]);
set(gca, 'Units', 'Pixels');
d = get(gca, 'Position');
im_width = round(d(3));
im_height = round(d(4));
t_height = 20;
t_width = 30;
br = 4;
l = 100;
y0 = im_height - t_height - 15;
for i = 1:N
    set(l1, 'XData', [t(i); t(i)]);
    set(l2, 'XData', [t(i); t(i)]);
    set(l3, 'XData', [t(i); t(i)]);
    set(l4, 'XData', [t(i); t(i)]);
    set(l5, 'XData', [t(i); t(i)]);
    theta = X(1, i)*pi/2;
    x0 = -X(3,i)*(im_width/2) + im_width/2;
    J = draw_pendulum(im_height, im_width, x0, y0, theta, br, l, t_height, t_width);
    imshow(J);
    drawnow
    frame = getframe(figh); % get the current frame
    im = frame2im(frame); % convert to normal image
    [A,map] = rgb2ind(im,256); %extract the image and color map
    if i == 1;
        imwrite(A,map,'animation.gif','gif', 'Loopcount',inf);
    else
        imwrite(A,map,'animation.gif','gif','WriteMode','append','DelayTime',1/24);
    end
end
rmpath('./pendulum');
rmpath('./util');