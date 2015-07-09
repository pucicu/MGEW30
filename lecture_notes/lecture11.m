%% Irregularly sampled time series analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct irregularly sampled time series

t = (-10000:2000)'; % time vector (e.g., Holocene)
N = length(t);

x1 = sin(2*pi*t/1000) + 0.5*randn(N,1);
x2 = cos(2*pi*t/1000) + 0.5*randn(N,1);
plot(t,x1,t,x2)

plot(x1,x2)






randn('seed',1);
t = (-10000:2000)'; % time vector (e.g., Holocene)
N = length(t);
x1 = zeros(N,1);
x2 = zeros(N,1);
rX1 = randn(N,1);
rX2 = randn(N,1);
x1(1) = rX1(1);
x2(1) = rX2(1);
a = 0.98;
k = .2;
for i = 2:N
   x1(i) = a * x1(i-1) + rX1(i);
   x2(i) = a * ((1-k) * x2(i-1) + k * x1(i-1)) + rX2(i) ;
end



% new (irregular) time axis
t1 = (-10000:5:1940)';
t2 = [-9800:8:-1500 -1480:4:2000]';
N1 = length(t1);
N2 = length(t2);
r1 = .7*randg(1,[N1 1]); 
r2 = .7*randg(1,[N2 1]); 
t1 = t1 + r1; t1 = sort(t1);
t2 = t2 + r2; t2 = sort(t2);

y1 = interp1(t,x1,t1); 
y2 = interp1(t,x2,t2); 

% now randomly remove some data points
rem = randi(length(t1), floor(length(t1)*0.75),1);
rem = [rem; find(isnan(y1))];
t1(rem) = []; y1(rem) = []; 

rem = randi(length(t2), floor(length(t2)*0.75),1);
rem = [rem; find(isnan(y2))];
t2(rem) = []; y2(rem) = [];

plot(t1,y1,'.-',t2,y2,'.-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolation method
% overlapping period
tStart = max(min(t1),min(t2));
tEnd = min(max(t1),max(t2));

% minimal and average sampling time
dtMean = mean([diff(t1); diff(t2)])
dtMin = min([diff(t1); diff(t2)])
dtMean = max([mean(diff(t1)); mean(diff(t2))])

% new regular time axis
tNew = tStart:dtMean:tEnd;

% interpolated time series
z1 = interp1(t1,y1,tNew);
z2 = interp1(t2,y2,tNew);

z1 = interp1(t1,y1,tNew,'spline');
z2 = interp1(t2,y2,tNew,'spline');

plot(t1,y1,'.:',t2,y2,'.:',tNew,z1,'.-',tNew,z2,'.-')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% kernel based interpolation
k = @(t,s)(exp(-(t.^2)./(2*s^2)) / sqrt(2*s^2 * pi));


s = .25 * mean(diff(t1));
z1k = zeros(length(tNew),1);
for i = 2:length(tNew)
    WL = k(t1 - tNew(i),s);
    z1k(i) = sum(WL .* y1) ./ sum(WL); 
end

plot(tNew,z1k,'.-',t1,y1,'.-')

plot(t,x1,'.-',tNew,z1k,'.-')
plot(tNew,z1,'.-',tNew,z1k,'.-')


s = 0.25 * mean(diff(t2));
z2k = zeros(length(tNew),1);
for i = 2:length(tNew)
    WL = k(t2 - tNew(i),s);
    z2k(i) = sum(WL .* y2) ./ sum(WL); 
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% smoothing with kernel

s = 2 * mean(diff(t1));
y1s = zeros(length(t1),1);
for i = 2:length(t1)
    WL = k(t1 - t1(i),s);
    y1s(i) = sum(WL .* y1) ./ sum(WL); 
end
plot(t1,y1s,t1,y1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% correlation

% original
[Corig lags_orig] = xcorr(x1,x2,1000,'coeff');
plot(lags_orig, Corig)

% interpolated
[Ci lags_i] = xcorr(z1,z2,round(1000/dtMean),'coeff');
plot(lags_i * dtMean, Ci)

% kernel interpolated
[Ck lags_k] = xcorr(z1k,z2k,round(1000/dtMean),'coeff');
plot(lags_k * dtMean, Ck)


% kernel based correlation

% normalize time series
y1n = zscore(y1);
y2n = zscore(y2);

% get all pairwise time distances 
t_dist = pdist2(t1, t2, @(x,y)(x-y)); t_dist = t_dist(:);

% product matrix of data
xy = pdist2(y1n,y2n, @(x,y)(x.*y)); xy = xy(:);

% apply lag
lags = (-1000:dtMean/2:1000);

C = zeros(length(lags),1);

s = 0.5 * dtMean;
for i = 1:length(lags)
    t_dist_lagged = t_dist + lags(i);
%    kernel_content = (t_dist_lagged <= 5*h) & (t_dist_lagged >= -5*h);
    % apply kernel
    WL = k(t_dist_lagged,s);
    Weight = sum(WL);
    % sum of the product of kernel with products x*y is the correlation
    C_ = sum(xy.*WL) ./ Weight;
    C(i) = C_;
end
    
plot(lags,C)



plot(lags_orig * 1, Corig,lags_i * dtMean, Ci,lags_k * dtMean, Ck,lags,C)

legend('orig','interp','kernel interp','kernel corr')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% autocorrelation


% normalize time series
y1n = zscore(y1);

% get all pairwise time distances 
t_dist = pdist2(t1, t1, @(x,y)(x-y)); t_dist = t_dist(:);

% product matrix of data
xy = pdist2(y1n,y1n, @(x,y)(x.*y)); xy = xy(:);

% apply lag
lags = (-10000:dtMean/1:10000);

C = zeros(length(lag),1);

s = 0.5 * dtMean;
for i = 1:length(lags)
    t_dist_lagged = t_dist + lags(i);
%    kernel_content = (t_dist_lagged <= 5*h) & (t_dist_lagged >= -5*h);
    % apply kernel
    WL = k(t_dist_lagged,s);
    Weight = sum(WL);
    % sum of the product of kernel with products x*y is the correlation
    C_ = sum(xy.*WL) ./ Weight;
    C(i) = C_;
end
    
plot(lags,C)


%% original
[Corig lags_orig] = xcorr(x1,x1,10000,'coeff');
plot(lags_orig, Corig)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% power spectrum
Fs = 1/mean(diff(lags));  % Sampling frequency

NFFT = 2^nextpow2(length(C));
Y = fft(C,NFFT)/length(C); % PSD
f = [Fs/2*linspace(0,1,NFFT/2+1)]'; % frequency vector

power = abs(Y);
power = (power(1:floor(NFFT/2)+1))'; % use only half 
power = power./sum(power);

plot(f,power)
xlim([0 .01])



%% original spectrum
[pxx fxx] = periodogram(x1,[],[],1);

plot(fxx,pxx)
xlim([0 .01])

%% spectrum from interpolated data
[pxx fxx] = periodogram(z1,[],[],1/dtMean);

plot(fxx,pxx)
xlim([0 .01])


%% Lomb Scargle

[Pxx Fxx Hxx] = lombscargle(t1,y1);


plot(Fxx,Pxx)
xlim([0 .01])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simple age modeling

clear

% load data
data = load('KRUM3_isotopes.dat');
dating = load('KRUM3_dates.dat');

% remove NaNs
data(isnan(data(:,4)),:) = [];
data(isnan(data(:,3)),:) = [];

% ages could be unsorted (correct to sorted ages)
dating(:,1)
dating = sortrows(dating,1);

% plot dates
plot(dating(:,3),dating(:,1),'.')
axis ij
ylabel('Depth (mm)'), xlabel('Age (yr BP)')

errorbar(dating(:,1),dating(:,3),dating(:,4))
xlabel('Depth (mm)'), ylabel('Age (yr BP)')

% interpolate a time axis for the original data (still irregularly sampled)
dataTime = interp1(dating(:,1),dating(:,3),data(:,1));

plot(dataTime, data(:,3))


% interpolate the data to a regularly sampled time axis

% create final time scale
fs = 4; % sampling time
t = dating(1,3):fs:dating(end,3); % equidistant time scale

% locations of an equally spaced time scale according to vector t
x = interp1(dating(:,3),dating(:,1),t,'lin');

% equally spaced data
o18 = interp1(data(:,1),data(:,4),x,'pchip');
c13 = interp1(data(:,1),data(:,3),x,'pchip');

% plot and compare time series of regular and irregular sampling
plot(t,o18)

plot(dataTime, data(:,4),'.-',t,o18,'.-')
xlabel('Age (yr BP)'), ylabel('d18O')
legend('irregularly sampled','regularly sampled')

%% Monte Carlo estimation of errors
M = 500; % number of realisations
o18 = zeros(length(t),M);
for i = 1:M
    dt = dating(:,4) .* randn(length(dating),1);
    x = interp1(dating(:,3) + dt,dating(:,1),t,'lin');
    o18(:,i) = interp1(data(:,1),data(:,4),x','pchip');
end

h = plot(t,o18);

% distribution of proxy values at a given time
idx = find(t >= 4000, 1); % index of time = 4000 yr BP
hist(o18(idx,:))

% plot median and 10% and 90% quantiles
o18median = median(o18,2);
o18q = quantile(o18',[.1 .9]); %'
plot(t, o18median, t, o18q,'r:')

legend('Median','10%, 90% quantile')
