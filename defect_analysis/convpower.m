function ret = convpower(x, n, rescale, cutoff)
if nargin < 3
    rescale = 1;
end
if nargin < 4
    cutoff = length(x);
end

prune = @(a) a(1:min(cutoff, length(a)));

ret = 0*x;
ret(1) = 1;
ret = prune(ret);

p = x;
while n > 0
    if rem(n,2) == 1
        ret = prune(conv(ret,p)) * rescale;
        n = n-1;
    end
    n = n/2;
    p = prune(conv(p,p)) * rescale;
end
end