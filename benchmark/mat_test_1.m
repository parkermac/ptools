% speed test

N = 1000;
M = 1000;

a = randn(N);


aa = sin(a);

tic

for ii = 1:M
    aa = sin(aa);
end

dt = toc;
disp(['took ', num2str(round(dt)), ' sec'])
