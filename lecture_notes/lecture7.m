%% Frequentist statistics 1

%% Empirical distributions
clear
randn('seed',1); % reset the random number generator to ensure reproducable random numbers
rand('seed',2); % reset the random number generator to ensure reproducable random numbers
x(:,1) = randn(10000,1); % normal distr. noise
x(:,2) = rand(10000,1); % uniform distr. noise
x(:,3) = randpink([10000,1],-2); % red noise
x(:,4) = randg(3,10000,1); % gamma dist. noise

for i = 1:4
   subplot(2,4,i)
   plot(x(:,i))
   subplot(2,4,i+4)
   hist(x(:,i),100)
end

%% Descriptors: central moments
% 1. moment = 0
% 2. moment = std dev
% 3. moment = skewness (after normalising)
% 4. moment = kurtosis (after normalising)
cm = inline('mean((x - mean(x)).^k)','x','k')

% central tendency
mean(x)
mode(x) % most frequent value
median(x) % median

% dispersion
std(x) % standard deviation
range(x) % range
diff(quantile(x,[.75 .25]))./sum(quantile(x,[.75 .25])) % quartile dispersion
var(x)./abs(mean(x)) % Fano factor

% distribution shape
skewness(x) % skewness
kurtosis(x) % kurtosis
kurtosis(x) - 3 % excess (kurt of norm distr = 3)

% quantiles, percentiles
quantile(x,[.25 .5 .75])
prctile(x,[25 50 75])

% summarising plot
boxplot(x)

% simple windowing with boxplot
clearrange
x = importdata('Mawmluh_village_precp2011.dat');
t = datenum(x.textdata(2:end,1),'dd.mm.yy');
x = x.data;
n = 14; % 14 days windows
y = round(linspace(1,length(x)/14,length(x)));

subplot(211)
plot(t,x), datetick
subplot(212)
boxplot(x,y*2)
xlabel('Week')


%% empirical distributions
% number of bins
clear
x = importdata('Mawmluh_village_precp2011.dat');
t = datenum(x.textdata(2:end,1),'dd.mm.yy');
x = x.data;
nBins = round(sqrt(length(x))) % number of bins

% bin centers
r = range(x);
cmin = min(x) + 0.5  * r/nBins;
cmax = max(x) - 0.5  * r/nBins;
bins = linspace(vmin, vmax, nBins);

% bin edges
emin = min(x);
emax = max(x) - r/nBins + eps;
binsLeft = linspace(emin, emax, nBins); % left bin edge
binsRight = linspace(emin+r/nBins, emax+r/nBins, nBins); % left bin edge


% count frequencies for each bin
h = zeros(nBins,1);
for i = 1:nBins
   h(i) = sum(x >= binsLeft(i) & x < binsRight(i));
end

% Matlab version
h2 = hist(x,nBins)

% new Matlab functions (since R2014b)
histogram(x)
h = histogram(x)
h.NumBins = 20
h.BinWidth = 50


h = histcounts(x)



%% Distributions
clear

% uniform
x = 0:10;
plot(x,unidpdf(x,5))

% alternatively pdf command
plot(x,pdf('unid',x,5))

% normal distr.
x = linspace(-5,20,1000);
m = 12; s = 1.16
p = normpdf(x,m,s);
plot(x,p)

p = p/sum(p);

% probability that x falls in the interval [m-2s m+2s]
fill([m-2*s m+2*s m+2*s m-2*s m-2*s],[0 0 max(p) max(p) 0],[.8 .8 .6])
hold on
plot(x,p)


idx = find((x >= m-2*s) .* (x <= m+2*s));
sum(p(idx))

% inverse pdf
norminv((1-0.95)/2,m,s)
norminv(1-(1-0.95)/2,m,s)

m+2*s,m-2*s

p = normspec([m-2*s m+2*s],m,s)
p = normspec([m-2*s m+2*s],m,s,'outside')



%% distribution fits

x = 2+3*randn(1000,1);
hist(x)
histfit(x)

[m s] = normfit(x)

%%
clear
x = importdata('Mawmluh_village_precp2011.dat');
t = datenum(x.textdata(2:end,1),'dd.mm.yy');
x = x.data;

[h bins] = hist(x);
bar(bins,h)

%%
m = expfit(x)

t = linspace(0,100,100);
plot(t,exppdf(t,m))

%%
p = gamfit(x)
plot(t,gampdf(t,p(1),p(2)))

%% Matlab functions:
p = mle(x,'dist','gamma') 

p = fitdist(x,'gamma') % new Matlab

gmdistribution
fitgmdist

%% Hypothesis tests

% Example: 2 coupled AR1 processes
clear
randn('seed',1); % reset the random number generator to ensure reproducable random numbers
N = 10000; % length of the time series
M = 1000; % number of realisations
x = zeros(N,M); y = zeros(N,M); % initialize vector with N times the value of zero
a = 0.98; % AR coefficient (for both systems the same)
k = .01; % coupling constant

% Calculate realisations of the coupled system
rX = randn(N,M); rY = randn(N,M); % create random numbers from Gaussian distribution
x(1,:) = rX(1,:); y(1,:) = rY(1,:); % initialze the first value in the time series with a random number
for i = 2:N
   x(i,:) = a * x(i-1,:) + rX(i,:);
   y(i,:) = a * ((1-k) * y(i-1,:) + k * x(i-1,:)) + rY(i,:) ;
end

% show time series
plot(x(:,1),y(:,1),'.')
plot(x(:,2),y(:,1),'.r')

% correlation between all combinations of time series
C_coup = diag(corr(x,y));

% Calculate realisations of the uncoupled system
k = 0;
for i = 2:N
   x(i,:) = a * x(i-1,:) + rX(i,:);
   y(i,:) = a * ((1-k) * y(i-1,:) + k * x(i-1,:)) + rY(i,:) ;
end
% correlation between all combinations of time series
C_uncoup = diag(corr(x,y));

% show results as histogram
clf
hist([C_coup C_uncoup],20) % only the coupled time series


% ROC
c_crit = linspace(-1,1,1000); % range of critical C-value
alpha = zeros(length(c_crit),1); % false positives
sens = zeros(length(c_crit),1); % true positives
for i = 1:length(c_crit)
   alpha(length(c_crit)-i+1) = sum(C_uncoup > c_crit(i))/length(C_uncoup);
   sens(length(c_crit)-i+1) = sum(C_coup > c_crit(i))/length(C_coup);
end
plot(alpha,sens)
xlabel('1-specificity')
ylabel('sensitivity')

% AUC
trapz(alpha,sens)

%% symmetric vs. asymmetric covariation
% example

% import data
del = load('delhi.dat');
jai = load('jaipur.dat');
plot(del(:,1),del(:,2),jai(:,1),jai(:,2))

subplot(211)
hist(del(:,2),linspace(0,100,20))
subplot(212)
hist(jai(:,2),linspace(0,100,20))

% corr
corr(del(:,2),jai(:,2))
corr(del(:,2),jai(:,2),'type','spearman')

% remove zeros
idx = find((del(:,2) == 0) + (jai(:,2) == 0) );
del(idx,:) = [];
jai(idx,:) = [];

subplot(211)
hist(del(:,2),linspace(0,100,20))
subplot(212)
hist(jai(:,2),linspace(0,100,20))

% corr
corr(del(:,2),jai(:,2))
corr(del(:,2),jai(:,2),'type','spearman')


% lognormal transformation
delLog = log(del(:,2));
jaiLog = log(jai(:,2));

plot(del(:,1),delLog,jai(:,1),jaiLog)

subplot(211)
hist(delLog,linspace(-5,5,20))
subplot(212)
hist(jaiLog,linspace(-5,5,20))

corr(delLog,jaiLog)


% differences
delDiff = diff(del(:,2));
jaiDiff = diff(jai(:,2));

plot(del(1:end-1,1),delDiff,jai(1:end-1,1),jaiDiff)

subplot(211)
hist(delDiff,linspace(-150,150,20))
subplot(212)
hist(jaiDiff,linspace(-150,150,20))

corr(delDiff,jaiDiff)


%% test distributions and hypothesis tests

%% t-test
% compare means of two distributions (normally distr.)
% (the difference between means follows t-distribution)
clear
N = 100;
% equal mean
randn('seed',1)
x1 = randn(N,1);
x2 = randn(N,1);

hist([x1 x2])

m1 = mean(x1); m2 = mean(x2);
s1 = var(x1); s2 = var(x2);

S = ((N-1) * s1 + (N-1) * s2) / (2*N-2)

T = sqrt(N^2/(2*N)) * abs(m1-m2)/S

% compare with t-distr.
clf
plot(linspace(-5,5,100),tpdf(linspace(-5,5,100),2*N-2))

% find crit. t-value
t_crit = tinv(1-0.05/2,2*N-2)

hold on
plot([t_crit t_crit],[0 .5],'r:') % critical t-value
plot([T T],[0 .5]) % measured t-value

% T is smaller than t_crit, i.e., HO cannot be rejected

% with Matlab function ttest2
[h p ci] = ttest2(x1,x2,.05)
mean(x1)-mean(x2)


% different mean
x1 = randn(N,1);
x2 = .5+randn(N,1);

hist([x1 x2])

[h p ci] = ttest2(x1,x2,.05)
[h p ci s] = ttest2(x1,x2,.05)
mean(x1)-mean(x2)

% t-test works only for normally distr. data!
% alternative: parameter free tests, e.g. Wilcoxon U-test

%% Wilcoxon U-test (= Mann-Whitney test)
clear
N = 100;
% equal mean
randn('seed',1)
x1 = 2*randg(2,[N,1]);
x2 = 7-randg(4,[N,1]);

hist([x1 x2])
mean(x1), mean(x2)

[h p ci] = ttest2(x1,x2,.05)
[p h s] = ranksum(x1,x2)


%% t-test for correlation
clear
N = 100;
% equal mean
randn('seed',1)
x1 = randn(N,1); x2 = randn(N,1); % uncoupled
C_uncoup = cor(x1,x2)

x1 = randn(N,1); x2 = 0.5 * randn(N,1) + 0.5 * x1; % coupled
C_coup = cor(x1,x2)


% t-value
T_coup = C_coup * sqrt(N-2) / sqrt(1-C_coup^2)
T_uncoup = C_uncoup * sqrt(N-2) / sqrt(1-C_uncoup^2)


% find crit. t-value
t_crit = tinv(1-0.10/2,2*N-2)




%% F-test
% compare variances

% same variance
x1 = randn(N,1);
x2 = 2+2*randn(N,1);

hist([x1 x2])

% F-value
s1 = var(x1);
s2 = var(x2); % s2 larger than s1

F = s2/s1

% crit F-value
finv(1-0.05/2, N-1, N-1)

% Matlab version
[h p ci] = vartest2(x2,x1,0.05)
[h p ci s] = vartest2(x2,x1,0.05)

% different variance
x1 = randn(N,1);
x2 = 2*randn(N,1);

hist([x1 x2])

[h p ci] = vartest2(x1,x2,0.05)



%% compare distributions
clear
x = importdata('Mawmluh_village_precp2011.dat');
x = x.data;
bins = linspace(0.1,300,10);
hx = hist(x,bins)

% fit exp distr
p = expfit(x)
hexp = exppdf(bins,p);

% rescale to fit to the data
hexp = hexp/sum(hexp) * sum(hx)

subplot(211)
bar(bins,hx)
subplot(212)
bar(bins,hexp)

% test if distr of x deviates from that of exp distr
X2 = sum((hx - hexp).^2 ./ hexp)

% crit. X2 value
chi2inv(1-0.05/2,10-(2+1))

%% KS-test
[cx xval] = ecdf(x); % emp. cum. distr of x
cexp = expcdf(xval,p); % exp. cum. distr

p = gamfit(x)
cexp = gamcdf(xval,p(1),p(2)); % exp. cum. distr

plot(xval,cx, xval, cexp)

KS = max(abs(cx-cexp))

% crit. KS-value
1.36/sqrt(length(x))

% Matlab function
kstest2

