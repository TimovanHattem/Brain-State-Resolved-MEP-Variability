function ccc = computeCCC(x, y)
    x = x(:); y = y(:);
    mx = mean(x);
    my = mean(y);
    vx = var(x, 1);   % population var
    vy = var(y, 1);
    cov_xy = mean((x - mx) .* (y - my));
    ccc = (2 * cov_xy) / (vx + vy + (mx - my)^2);
end