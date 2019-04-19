%% Simulation 1: comparing meansurements of different measurements
% In this simulation, we evaluate the perfomance of GI and range of ''CI''
% when we are comparing two different measurements with different scale. We
% denote the two measurements by X and Y. X and |Y| follow a beta
% distribution parametrized by alpha and beta, where alpha =1 and beta
% follows from a gamma distribution. The sign of Y is randomly chosen with
% equal probability to be + or -. Consequently, the support of X is [0,1]
% while the support of Y is [-1,1].

rng(0)
E = 10; %number of experiemtns
n = 1000; %number of subjects in one experiment;
alpha = 1; % first parameter 
betax = exprnd(1,E,1); % second parameter for x
betay = 4*ones(E,1);
X = zeros(E,n);
Y = zeros(E,n);
UpperCI = zeros(E,2);
lowerCI = zeros(E,2);

for e=1:E
    X(e,:) = betarnd(alpha,betax(e),1,n);
    Y(e,:) = betarnd(alpha,betay(e),1,n)*sign(normrnd(1,n));
    UpperCI(e,1) = quantile(X(e,:),0.975);
    lowerCI(e,1) = quantile(X(e,:),0.025);
    UpperCI(e,2) = quantile(Y(e,:),0.975);
    lowerCI(e,2) = quantile(Y(e,:),0.025);
end
%% Calculate GI
l = max(UpperCI(:,1)) - min(lowerCI(:,1));
O = min(UpperCI(:,1)) - max(lowerCI(:,1));
GIx  = (l-O)/sqrt(E);
GIx =exp(-5*GIx);

l = max(UpperCI(:,2)) - min(lowerCI(:,2));
O = min(UpperCI(:,2)) - max(lowerCI(:,2));
GIy  = (l-O)/(2*sqrt(E));
GIy =exp(-5*GIy);
%% Calculate ''CI range''
CIrangeX = quantile(X(:),0.975)-quantile(X(:),0.025);
CIrangeY = quantile(Y(:),0.975)-quantile(Y(:),0.025);
%%
g = repmat(1:E,n,1);
subplot(1,2,1)
boxplot(X(:),g(:),'Notch','off','Whisker',0,'Symbol','');
%bar(lowerCI(:,1),UpperCI(:,1));
ylim([-1,2]);
subplot(1,2,2)
boxplot(Y(:),g(:),'Notch','off','Whisker',0,'Symbol','');
%bar(lowerCI(:,2),UpperCI(:,2));
ylim([-1,1]);

