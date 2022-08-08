function F = fullGaussian(x,xdata)
% [Ib, I0, x0, sx, y0, sy, theta]

a = cos(x(7)).^2/2/x(4).^2+sin(x(7)).^2/2/x(6).^2;
b = -sin(2*x(7))/4/x(4).^2+sin(2*x(7))/4/x(6).^2;
c = sin(x(7)).^2/2/x(4).^2+cos(x(7)).^2/2/x(6).^2;

F = x(1) + x(2)*exp(-(...
    a*( xdata(:,:,1) - x(3) ).^2 + ...
    2*b*( xdata(:,:,1) - x(3) ) .* ( xdata(:,:,2) - x(5) ) + ...
    c*( xdata(:,:,2) - x(5) ).^2 ));

