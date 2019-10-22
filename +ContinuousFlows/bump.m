function y = bump(x, interval, delta)
% BUMP Bump function zero outside interval [a,b], and = 1 within
%      [a+delta,b-delta].
%%     function y = bump(x, interval, delta)


a = interval(1);
b = interval(2);

y = nan(size(x));
y( x < a | x > b) = 0;
y( x > (a+delta) & x < (b-delta) ) = 1;

lsel = x >= a & x <= a+delta;
rsel = x >= b-delta & x <= b;

y( lsel ) = (x(lsel) - a)/delta;
y( rsel ) = (b - x(rsel))/delta;
