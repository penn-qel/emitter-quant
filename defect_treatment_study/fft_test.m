n = 360;
d = 0.97;
x = linspace(0, 2*pi, n);
y = cos(2*x + d) + wgn(1,n,-10) + cos(0.357*x.^2 + 0.789);
z = fft(y);
z3 = [z(1:3) zeros(1,19) z(23:end)];
%plot(x, y - ifft(z3))

n2 = 1000*n;
x2 = linspace(0, 2*pi, n2);
y2 = interp1(x, y, x2, 'spline');
z2 = fft(y2);
d-angle(z(3))
d-angle(z2(3))