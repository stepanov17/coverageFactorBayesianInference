clear;
clc;


function [maxR, appr, d] = sAppr(x, y, bp, plotFlag, plotCol)

    nx = length(x);
    xStart = x(1); xEnd   = x(nx);

    [appr, d] = lsq_splin(x, y, bp);

    kRound = 1.e3;
    appr = (1 ./ kRound) * round(kRound * appr);
    d = (1 ./ kRound) * round(kRound * d);

    yAppr = interp(x, bp, appr, d, "natural");

    if (plotFlag > 0) then
        plot(x, y, plotCol + '.')
        plot(x, yAppr, plotCol)
        // 'x'
        yAppr2 = interp(bp, bp, appr, d, "natural");
        plot(bp, yAppr2, plotCol + 'x')
    end

    r = abs(y - yAppr) ./ y;
    maxR = max(r * 100); // max rel. error, %

endfunction


// conv(N, E(2)), p0 = 99, n = 4
x = [0.    0.1   0.2   0.3   0.4   0.5   0.75  1.    1.25  1.5   1.75  2.    2.25  2.5   2.75  3.    3.25  3.5   3.75  4.    4.5   5.    5.5   6.    7.    8.    9.    10.  ];
y = [3.372 3.371 3.362 3.342 3.319 3.287 3.194 3.087 2.979 2.887 2.816 2.759 2.717 2.687 2.668 2.650 2.638 2.627 2.623 2.613 2.605 2.601 2.598 2.594 2.589 2.584 2.582 2.582];

xMax = x(length(x));
h = 0.005;
nh = xMax / h;

iOpt = [];
R = 1.e8;

nIter = 1;

n1 = 0.35 * nh; // nh - 1;
n2 = 0.70 * nh; // nh - 2;

// quasi-optimal breakpoints, argmin(maxR)
for i1 = 1 : n1
    for i2 = i1 + 1 : n2
        breakpoints = [0., i1 * h, i2 * h, xMax];
        [maxR, yi, di] = sAppr(x, y, breakpoints, -1, 'k');
        if (maxR < R) then
            R    = maxR;
            iOpt = [i1 i2];
        end
        printf(">> %d\t[%d %d]\t%f\n", nIter, i1, i2, R);
        nIter = nIter + 1;
    end
end

breakpoints = [0, h * iOpt, xMax]
[maxR, yi, di] = sAppr(x, y, breakpoints, 1, 'k');
maxR

splineData = [breakpoints; yi; di]
printf("maxR = %.2f", maxR)
