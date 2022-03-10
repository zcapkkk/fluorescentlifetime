function loss = expdec3eval(x, data)
    t = data(:,1);
    y = data(:,2);
%     y = y/max(y);
    fit = x(1)*exp(-t/x(4))...
        +x(2)*exp(-t/x(5))...
        +x(3)*exp(-t/x(6))+x(7);
    loss = sum((y - fit).^2./y)/(length(t)-length(x));
%     loss = sum((y - fit).^2./y)/(length(t)-length(x));
%     loss = sum((y - fit).^2);
end
