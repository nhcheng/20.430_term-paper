%% Figure 1
clc
clear
syms x
txt = [];
lambda_hist = linspace(4,28,7);
for i = 1:length(lambda_hist)
lambda = lambda_hist(i);
analy_sp(lambda,length(lambda_hist),i)
legendInfo{i}=[strcat('\lambda = ',' ', num2str(lambda))];
end
legend(legendInfo)

%% Figure 2
clc
clear
syms x
txt = [];
lambda_hist = linspace(4,28,7);
for i = 1:length(lambda_hist)
lambda = lambda_hist(i);
analy(lambda,length(lambda_hist),i)
legendInfo{i}=[strcat('\lambda = ',' ', num2str(lambda))];
end
legend(legendInfo)

%% Figure 3
clc
clear
close all
c_hist = logspace(-1,3,9);
for i = 1:length(c_hist)
c = c_hist(i);
flow(c,length(c_hist),i)
hold on
end
xline(0.7,'--k','LineWidth',2.5,'DisplayName','20% MFI');

%% Figure 4
clc
clear
close all
lam_hist = linspace(4,28,7); 
lam_hist = [0.004,0.01,0.03,0.09,0.2,0.6,1]
for i = 1:length(lam_hist)
c = 2*fliplr(lam_hist(i))^(-1/2);
flow_c(c,length(lam_hist),i,lam_hist)
hold on
end
xline(0.7,'--k','LineWidth',2.5,'DisplayName','20% MFI');

%% Figure 5
syms x
y = 2*x^(-1/2)
fplot(y,'LineWidth',2)
hold on
y = 3.8*x^(-1/2)
fplot(y,'LineWidth',2)
set(gca, 'xscale','log')
set(gca, 'yscale','log')
xlim([1e-2 1])
ylim([1e-1 1000])
grid on 
hold on
xlabel('IL-2R\alpha^+ density per 10^3\mum^3')
ylabel('\lambda_{diameter}')
legend('Simplified model','Oyler-Yaniv et al.')
title('Scaling relationship of simplified model and Oyler-Yaniv et al.')

%%
function flow(c,l,i)
color = winter(l);
syms x
fun = 4/(1+6.5/x)+0.1;
x = c*ones(1,100000);
ffun = matlabFunction(fun);
res = ffun(x) +0.1; % 0.1 as background
res = abs(randn(1,length(res)).*res);
[N,edges] = histcounts(log10(res));
[N,edges]  = histcounts(res,10.^edges);
plot(conv(edges, [0 1], 'valid'), 100*N/max(N),'DisplayName',strcat(['\alpha = ',' '...
    ,num2str(log10(c))]),'LineWidth',2.5,'color',color(i,:))
set(gca, 'xscale','log')
xlim([0.01 50])
xlabel('pSTAT5 (a.u.)')
ylabel('% of Max')
legend('Location','northeastoutside')
grid on
title('Model of pSTAT5 response of CD4^+IL-2R\alpha^+ consuming cells to [IL-2]=10^\alpha')
end

function flow_c(lambda,l,i,lam_hist)
% 2mm length = 2000 micrometer = 200 cell diameters
% Cell: 10 micrometer
color = autumn(l);
syms x
%y = linspace(1,200,1e8);
y = linspace(1,200,10000000);
c = 200*exp(-1/lambda*x); %200pm
%c = 400*(lambda/x)^2
f = 4/(1+6.5/c);
ff = matlabFunction(f);
res = (ff(y))+0.1; %background
res = abs(randn(1,length(res)).*res);
[N,edges] = histcounts(log10(res));
[N,edges]  = histcounts(res,10.^edges);
plot(conv(edges, [0 1], 'valid'),100*N/max(N),'DisplayName',strcat(['[IL-2R\alpha^+] = ',' '...
    ,num2str(lam_hist(i))]),'LineWidth',2.5,'color',color(i,:))
set(gca, 'xscale','log')
xlim([0.01 20])
xlabel('pSTAT5 (a.u.)')
ylabel('% of Max')
legend('Location','northeastoutside')
grid on
title('Distribution of pSTAT5 in CD4^+IL-2R\alpha^+ in clusterwell with varying IL-2R\alpha^+ density')
end

function analy(lambda,l,i)
color = flipud(summer(l));
syms x
f = 4/(1+6.5/(200*exp(-1/lambda*x)));
g = finverse(f);
h = -diff(g)/200
%h = subs(h,x,1/4*x)
figure(1)
fplot(h,'LineWidth',2.5,'color',color(i,:))
% y = rand(10000,1)*200;
% c = 100*exp(-0.98*x);
grid on 
hold on
xlim([-0.01 4/(1+6.5/200)])
ylim([0 1])
ylabel('Probability density')
xlabel('Mean fluorescence intensity (MFI)')
title('Probability distribution of pSTAT5 MFI in CD4^+IL-2R\alpha^+ in clusterwell with varying \lambda_{diameter}')
legend('Location','northeastoutside')
end

function analy_sp(lambda,l,i)
color = flipud(summer(l));
syms x
f = 4/(1+6.5/(200*(exp(-1/lambda*x)-exp(-1/lambda*200))));
g = finverse(f);
h = -diff(g)/200;
%h = subs(h,x,1/4*x)
figure(1)
fplot(f,'LineWidth',2.5,'color',color(i,:))
% y = rand(10000,1)*200;
% c = 100*exp(-0.98*x);
ylim([-0.01 4/(1+6.5/200)])
ylim([1e-3 4/(1+6.5/200)])
xlim([0 200])
grid on 
hold on
xlabel('Distance (Cell diameters)')
ylabel('Mean fluorescence intensity (MFI)')
title('Spatial variation of pSTAT5 MFI in CD4^+IL-2R\alpha^+ in clusterwell with varying \lambda_{diameter}')
legend('Location','northeastoutside')
end