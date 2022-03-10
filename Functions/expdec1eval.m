function loss = expdec1eval(x, data)
    t = data(:,1);
    y = data(:,2);
%     y = y/max(y);
    fit = x(1)*exp(-t/x(2))+x(3);
    loss = sum((y - fit).^2./y)/(length(t)-length(x));
    %loss = sum((y - fit).^2);
end
