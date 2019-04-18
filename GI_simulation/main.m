%% Simulation 1: comparing meansurements of different measurements
% In this simulation, we evaluate the perfomance of GI and range of ''CI''
% when we are comparing two different measurements with different scale. We
% denote the two measurements by X and Y. X and |Y| follow a beta
% distribution parametrized by alpha and beta, where alpha =1 and beta
% follows from a gamma distribution. The sign of Y is randomly chosen with
% equal probability to be + or -. Consequently, the support of X is [0,1]
% while the support of Y is [-1,1].


E = 6; %number of experiemtns
n = 100; %number of subjects in one experiment;
alpha = 1; % first parameter 
betax = exprnd(10,E,1); % second parameter for x
betay = exprnd(5,E,1);
X = zeros(E,n);
Y = zeros(E,n);
UpperCI = zeros(E,2);
lowerCI = zsros(E,2);

for e=1:E
    X(e,:) = betarnd(alpha,betax(e),1,n);
    Y(e,:) = betarnd(alpha,betay(e),1,n)*sign(normrnd(1,n));
end
g = repmat(1:E,n,1);
%%
boxplot(Y(:),g(:),'Notch','off');
ylim([-2,2]);







