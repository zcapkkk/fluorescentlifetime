function y = expdec(params, t)
n = length(params);
y = zeros(length(t),1);
for i = 1:2:n-1
    y = y + params(i)*exp(-t/params(i+1));
end
y=y+params(n);
