function y = expdec(params, t)
    n = length(params);
    splithalf = (n-1)/2;
    y = zeros(length(t),1);
    amplitudes = params(1:splithalf);
    periods = params(splithalf+1:end-1);
    for i = 1:splithalf
        y = y + amplitudes(i)*exp(-t/periods(i));
    end
    y=y+params(end);
end
