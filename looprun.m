addpath('./Functions');
close all;

amplitudes0 = [0.8, 0.1, 0.1];
lifetimes0 = [0.4, 1.5, 4.4];

times = (1:5:60)*60;
monochi = zeros(1,length(times));
bichi = zeros(1,length(times));
trichi = zeros(1,length(times));

for i = 1:length(times)
   [monodata, bidata, tridata] = expeval(amplitudes0, lifetimes0, times(i));
   monochi(i) = monodata(end-1);
   bichi(i) = bidata(end-1);
   trichi(i) = tridata(end-1);
%    disp(monodata(end));
%    disp(bidata(end));
%    disp(tridata(end));
end

figure;

times = times/60;
plot(times, monochi,'--x');
hold on;
plot(times, bichi,'--x');
plot(times, trichi,'--x');
hold off;
legend("Mono","Bi","Tri");
grid on;
title("Reduced Chi Square for different acquisition times");
xlabel("Acquisition times (mins)");
ylabel("Reduced Chi Square");