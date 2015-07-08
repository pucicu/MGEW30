%% hypothesis testing with resampling, shuffling, block shuffling
% constructing empirical test distr.
clear
x = load('weather_potsdam.dat');
% remove NaNs
x(isnan(x(:,2)),:) = [];
plot(x(:,1),x(:,2))
datetick('keeplimits')

y = 1890:10:2015;
a = zeros(length(y)-1,3);
for i = 1:length(y)-1
    idx = find((x(:,1) >= datenum(y(i),1,1)) & (x(:,1) < datenum(y(i+1),1,1)));
    a(i,:) = arfit(x(idx,2),2); 
end
nw = length(idx); % length of window

plot(y(1:end-1),a(:,1))

%% random shuffling (bootstrap method)
M = 100;
as = zeros(M,3);
h = waitbar(0,'Surrogate test');
for j = 1:M, waitbar(j/M)
    xs = x(randperm(length(x),nw),2);
    as(j,:) = arfit(xs,2); 
end
close(h)

clf
hist(as(:,1))

% quantiles as CI
ci = quantile(as(:,1),[.05 .95])

% plot results together with CI
clf
plot(y(1:end-1),a(:,1))
hold on
plot(y(1:end-1), ones(length(y)-1,1)*ci(1),'r:',y(1:end-1), ones(length(y)-1,1)*ci(2),'r:')

%% block bootstrap
M = 100;
as = zeros(M,3);
h = waitbar(0,'Surrogate test');
w = 364/2; % length of block
for j = 1:M, waitbar(j/M)
    idx = randperm(length(x)-w,nw)';
    idx = repmat(idx,1,w) + repmat(0:w-1,length(idx),1);
    idx = idx';
    idx = idx(:);
    xs = x(idx(1:nw),2);
    as(j,:) = arfit(xs,2); 
end
close(h)

hist(as(:,1))

% quantiles as CI
ci = quantile(as(:,1),[.05 .95])

% plot results together with CI
clf
plot(y(1:end-1),a(:,1))
hold on
plot(y(1:end-1), ones(length(y)-1,1)*ci(1),'r:',y(1:end-1), ones(length(y)-1,1)*ci(2),'r:')



%% phase shuffling (takes too long and is actually not appropriate for this application)
M = 100;
as = zeros(length(y)-1,3,M);
h = waitbar(0,'Surrogate test');
xfft = fft(x(:,2));
for j = 1:M, waitbar(j/M)
    idx = randperm(length(x));
    xffts = complex(real(xfft), imag(xfft(idx)));
    xs = ifft(xffts,'symmetric');
    for i = 1:length(y)-1
        idx = find((x(:,1) >= datenum(y(i),1,1)) & (x(:,1) < datenum(y(i+1),1,1)));
        as(i,:,j) = arfit(xs(idx),2); 
    end
end
close(h)

clf
hist(squeeze(as(:,1,:)))

% quantiles as CI
ci = zeros(length(y)-1,2);
for i = 1:length(y)-1
   ci(i,:) = quantile(squeeze(as(i,1,:)),[.05 .95])
end


% plot results together with CI
clf
plot(y(1:end-1),a(:,1))
hold on
plot(y(1:end-1), ci(:,1),'r:',y(1:end-1), ci(:,2),'r:')


%% Matlab bootstr
fano = @(x)[var(x)./abs(mean(x))]; % Fano factor

a = zeros(length(y)-1,1);
for i = 1:length(y)-1
    idx = find((x(:,1) >= datenum(y(i),1,1)) & (x(:,1) < datenum(y(i+1),1,1)));
    a(i) = fano(x(idx,2)); 
end
nw = length(idx); % length of window

as = bootstrp(nw,@(x)[var(x)./abs(mean(x))],x(:,2))

% quantiles as CI
ci = quantile(as(:,1),[.05 .95])

% plot results together with CI
clf
plot(y(1:end-1),a)
hold on
plot(y(1:end-1), ones(length(y)-1,1)*ci(1),'r:',y(1:end-1), ones(length(y)-1,1)*ci(2),'r:')


%% Multiple test effects, correction schemes

% test the same data set, but many variables
a = zeros(length(y)-1,4);
for i = 1:length(y)-1
    idx = find((x(:,1) >= datenum(y(i),1,1)) & (x(:,1) < datenum(y(i+1),1,1)));
    a(i,1) = fano(x(idx,2)); % temp
    a(i,2) = fano(x(idx,6)); % prcp
    a(i,3) = fano(x(idx,7)); % rH
end

% apply bootstrap procedure on 3 variables
as = bootstrp(nw,@(x)[var(x)./abs(mean(x))],x(:,2));
as(:,2) = bootstrp(nw,@(x)[var(x)./abs(mean(x))],x(:,6));
as(:,3) = bootstrp(nw,@(x)[var(x)./abs(mean(x))],x(:,7));



% plot results together with CI (consider only exceeding dispersion, i.e., one-sided test)
alpha = 0.05;
ci = quantile(as,1-alpha)
for i = 1:3
    subplot(3,1,i)
    plot(y(1:end-1),a(:,i))
    hold on
    plot(y(1:end-1), ones(length(y)-1,1)*ci(1,i),'r:',y(1:end-1), ones(length(y)-1,1)*ci(2,i),'r:')
end

% we have 12x3 = 36 tests
size(a)
m = prod(size(a))

% estimate p-value
p = zeros(length(a),3);
for i = 1:3
    [ec xval] = ecdf(as(:,i)); % create cum. empir. test distr.
    for j = 1:length(a)
       idx = find(xval < a(j,i),1,'last'); % find upper
       if isempty(idx)
           p(j,i) = 1;
       else
           p(j,i) = 1-ec(idx);
       end
    end
end

%% Benjamini-Hochberg correction
p = p(:);
[ps idx] = sort(p);

% number of th H0 to be rejected
k = sum(ps < (1:m)' * alpha/m)

% prepare mask matrix for sign. results
b = zeros(size(a));
b(idx(1:k)) = 1;

% mark the "really" signifcant values by a "*" 
for i = 1:3
    subplot(3,1,i)
    plot(y(find(b(:,i))),a(find(b(:,i)),i),'*')
end

