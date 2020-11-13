clear;
clc;


function [maxR, appr, d] = sAppr(x, y, bp)

    nx = length(x);
    xStart = x(1); xEnd   = x(nx);

    [appr, d] = lsq_splin(x, y, bp);

    kRound = 1.e3;
    appr = (1 ./ kRound) * round(kRound * appr);
    d = (1 ./ kRound) * round(kRound * d);

    yAppr = interp(x, bp, appr, d, "natural");

    plot(x, y, 'k.')
    plot(x, yAppr, 'r')

    // 'x'
    yAppr2 = interp(bp, bp, appr, d, "natural");
    plot(bp, yAppr2, 'rx')

    r = abs(y - yAppr) ./ y;
    maxR = max(r * 100); // percent

endfunction


DATA = fscanfMat("N_95_5.txt");
sz = size(DATA);

nData = sz(1);
nx = sz(2);

x = DATA(2, :);
x = x(2 : nx);

breakpoints = [0., 0.75, 1.75, 3.5, 10.];

R = [];
splineData = [];

for i = 3 : nData
    y = DATA(i, :);
    y = y(2 : nx);
    [maxR, yi, di] = sAppr(x, y, breakpoints);
    splineData = [splineData; yi; di];
    R = [R, maxR];
end

printf("max rel. err = %3.2f %%\n", max(R))

splineData
