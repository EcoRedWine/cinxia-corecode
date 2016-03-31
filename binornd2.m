function r = binornd2(n,p)
if isscalar(p)
p = p * ones(size(n));
elseif size(n) ~= size(p)
    error('length of n and p different in binornd2')
end

r = zeros(size(n));

    for i = 1:max(n(:))
        k = find(n >= i);
        r(k) = r(k) + (rand(size(k)) < p(k));
    end
    r(p < 0 | 1 < p | n < 0 | round(n) ~= n) = NaN;
end