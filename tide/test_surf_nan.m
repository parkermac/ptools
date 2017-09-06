% code to test surf with an array containing nan's

x = linspace(0,10,11);
y = x;
[X, Y] = meshgrid(x,y);

Z = X.^2 + Y.^2;

Z(4,4) = NaN;

sh = surf(X, Y, Z);

shading interp


