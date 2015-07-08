%% lorenz attractor

% integrate Lorenz system
t2=1:.01:100;

[t,z]=ode45(@l,[t2],[2.5798 5.5037 24.2663]);

[t,z]=ode45('l',[t2],[2 2 2]);

z(1:4000,:) = []; % remove transients
t(1:4000) = [];
N = length(z); 

% time series
plot(t,z)

% phase space
plot3(z(:,1),z(:,2),z(:,3))

%% phase space reconstrucion using one component

% autocorrelation
clf
maxlag = 500;
C = xcorr(z(:,1),maxlag,'coeff');
plot(-maxlag:maxlag,C,'.-')

% correlation time
hold on
plot([-maxlag maxlag],[1/exp(1) 1/exp(1)],'r')

find(C(maxlag+1:end) < 1/exp(1),1)

% show reconstructed phase space
clf
tau = 10;
plot3(z(1:end-2*tau,1),z((1+tau):end-tau,1),z((1+2*tau):end,1))

%% correlation integral & dimension
e=10.^[-2:.15:1.5]; m = 1:8;
h = waitbar(0,'Calculation');
for i = 1:length(e);
    for j = 1:length(m); waitbar((i*length(m)+j)/(length(e)*length(m)))
        X = crp(z(1:2:5000,1),m(j),1,e(i),'max','sil');
        RR(i,j) = mean(X(:));
    end
end
close(h)

% correlationsum
loglog(e,RR)
xlabel('Threshold \epsilon'), ylabel('C(\epsilon)')
legend(num2str(m')), grid

% locale slope
i = 3:10;
for j = 1:length(m)
    dummy = polyfit(log(e(i))',log(RR(i,j)),1);
    d(j) = dummy(1);
end

plot(m,d,'.-')
xlabel('Embedding dimension m'); ylabel('D_2')
grid


% mutual information
m = mi(z(:,1),'sil');

% false nearest neighbours
fnn(z(:,1))

