function F = Gaussian1D(x,xdata)
% [Ib, I0, z0, sz]

F = x(1) + x(2)*exp(-(xdata-x(3)).^2/2/x(4).^2 );

