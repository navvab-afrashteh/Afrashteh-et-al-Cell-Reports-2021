function Theta = angleUV(u,v)
% calculates the angle between 2D vectors u and v
x1 = real(u); y1 = imag(u);
x2 = real(v); y2 = imag(v);

Theta = atan2d(x1*y2-y1*x2,x1*x2+y1*y2);

