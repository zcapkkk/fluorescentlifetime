addpath('./Functions');
close all;
format long; 

% amplitudes0 = [0.66,0.33];
% lifetimes0 = [0.23,0.54];

% Demo
% amplitudes0 = [0.8, 0.1, 0.1];
% lifetimes0 = [0.4, 1.5, 4.4];

% mCherry
% lifetimes0 = [1.47, 2.77];
% amplitudes0 = [0.84, 0.16];

% Yellow Camelon
lifetimes0 = [0.67,2.2,3.57];
amplitudes0 = [0.45,0.36,0.19];


lifetimes0 = sort(lifetimes0);
amplitudes0 = sort(amplitudes0);


times = (1:120)*60;
monochi = zeros(1,length(times));
bichi = zeros(1,length(times));
trichi = zeros(1,length(times));
monoerr = zeros(3,length(times));
bierr = zeros(5,length(times));
trierr = zeros(7,length(times));
fdistvals = zeros(2,length(times));
monoparams = zeros(3,length(times));
biparams = zeros(5,length(times));
triparams = zeros(7,length(times));

for i = 1:length(times)
   [monodata, bidata, tridata, err1, err2, err3, ft12, ft23] = expeval(amplitudes0, lifetimes0, times(i));
   rc1 = monodata(end-1);
   rc2 = bidata(end-1); 
   rc3 = tridata(end-1);
   
   monochi(i) = rc1;
   bichi(i) = rc2;
   trichi(i) = rc3; 
   monoerr(:,i) = err1;
   bierr(:,i) = err2;
   trierr(:,i) = err3;
   % threshold is 0.05, 
   % if less reject H0 and accept more complex models
   fdistvals(:,i) = [ft12, ft23];
   monoparams(:,i) = monodata(1:3);
   biparams(:,i) = bidata(1:5);
   triparams(:,i) = tridata(1:7);
end

% data processing

% generate amplitude arrays
monoamps = monoparams(1,:);
biamps = biparams(1:2,:);
biamps = biamps./sum(biamps);
triamps = triparams(1:3,:);
triamps = triamps./sum(triamps);
% generate amplitude errors
monoampserr = monoerr(1,:);
biampserr = bierr(1:2,:)./biamps;
triampserr = trierr(1:3,:)./triamps;

% generate lifetime arrays
monolfts = monoparams(2,:);
bilfts = biparams(3:4,:);
trilfts = triparams(4:6,:);
% generate lifetime errors
monolftserr = monoerr(2,:);
bilftserr = bierr(3:4,:);
trilftserr = trierr(4:6,:);

times = times/60;

figure("Name","Reduced Chi Square Plot all");
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

figure("Name","Reduced Chi Square Plot bi/tri");
plot(times, bichi,'--x');
hold on;
plot(times, trichi,'--x');
hold off;
legend("Bi","Tri");
grid on;
title("Reduced Chi Square for different acquisition times");
xlabel("Acquisition times (mins)");
ylabel("Reduced Chi Square");


% figure('Name','Errors');
% subplot(211);
% plot(times,monolfts,'-.');
% hold on;
% for i = 1:2
%     plot(times, biamps(i,:), '-.');
% end
% for i = 1:3
%     plot(times, triamps(i,:), '-.');
% end
% hold off;
% grid on;
% title("Amplitude parameters");
% xlabel("Acquisition times (mins)");
% ylabel("Error in amplitude estimates");
% subplot(212);
% plot(times,monolfts,'-x');
% hold on;
% for i = 1:2
%     plot(times, bilfts(i,:), '-x');
% end
% for i = 1:3
%     plot(times, trilfts(i,:), '-x');
% end
% hold off;
% grid on;
% title("Lifetime parameters");
% xlabel("Acquisition times (mins)");
% ylabel("Error in lifetime estimates");

% f test plots
figure("Name","F test");
plot(times, fdistvals(1,:),'LineWidth',1.5);
hold on;
plot(times, fdistvals(2,:),'LineWidth',1.5);
hold off;
grid on;
legend("bi to mono","tri to bi");
title("F test results");
xlabel("acquisition time (mins)");



% comparison with actual
n = length(amplitudes0);


if n == 2
    ampdiff = sort(biamps,1);
    lifediff = sort(bilfts,1);
elseif n == 3
    ampdiff = sort(triamps,1);
    lifediff = sort(trilfts,1);
end

for i = 1:n
    ampdiff(i,:) = abs(ampdiff(i,:) - amplitudes0(i))/amplitudes0(i);
    lifediff(i,:) = abs(lifediff(i,:) - lifetimes0(i))/lifetimes0(i);
end

% figure('Name','Difference');
% subplot(211);
% hold on;
% for i = 1:n
%     plot(times, ampdiff(i),'-x');
% end
% title("Percentage error from actual amplitudes")
% xlabel("acquisition time (mins)");
% ylabel("absolute percentage");
% grid on;
% 
% subplot(212);
% for i = 1:n
%     plot(times, lifediff(i),'-x');
% end
% title("Percentage error from actual lifetimes");
% xlabel("acquisition time (mins)");
% ylabel("absolute percentage");
% grid on;

figure("Name","Deviation from actual");
plot(times, mean(lifediff,1)*100,'-.','LineWidth',2);
hold on;
for i = 1:n
    plot(times, lifediff(i,:)*100,'LineWidth',1.5);
end
yline(1,'--');
hold off;
if n == 2
    legend("average","t1","t2");
elseif n == 3
    legend("average","t1","t2","t3");
end
title("Percentage error from actual lifetime value");
xlabel("acquisition time (mins)");
ylabel("percentage (%)");
grid on;


figure("Name","Hessian lifetimes error");
plot(times, monolftserr,'LineWidth',2);
hold on;
for i = 1:2
    plot(times, bilftserr(i,:),'LineWidth',2);
end
for i = 1:3
    plot(times, trilftserr(i,:),'LineWidth',2);
end
hold off
grid on;
legend('Mono','Bi t1','Bi t2','Tri t1','Tri t2','Tri t3');
title("Hessian Lifetime Error");
xlabel("acquisition time (mins)");
ylabel("error (ns)");
    
