%% Correlation analysis of coupled AR processes
% In this example we analyze the correlation between two
% coupled AR(1) processes in dependence of the coupling
% strength.
%
% Homework 3

% (coding: Norbert Marwan, 5/2015)

%% Initialisation
% Some first settings.
clear
randn('seed',1); % reset the random number generator to ensure reproducable random numbers
N = 10000; % length of the time series
x = zeros(N,1); % initialize vector with N times the value of zero
y = zeros(N,1); % initialize vector with N times the value of zero
rX = randn(N,1); % create random numbers from Gaussian distribution
rY = randn(N,1); % create random numbers from Gaussian distribution
x(1) = rX(1); % initialze the first value in the time series with a random number
y(1) = rY(1); % initialze the first value in the time series with a random number

%% Parameter
a = 0.98; % AR coefficient (for both systems the same)
L = 500; % number of different coupling constants to be analysed
maxLag = 100; % maximum lag for cross correlation

% correlation coefficient
C = zeros(2*maxLag+1,L);
cnt = 1; % counter

%% Calculation of correlation coefficient for differen coupling
h = waitbar(0, 'Calculate cross correlation'); % create a waitbar
for k = linspace(0,1,L); % coupling constant
    waitbar(cnt/L) % update waitbar
    % Calculate the coupled system
    % The system _x_ is a standard AR(1), the system _y_ is coupled with
    % system _x_ by a coupling constant _k_. In order to preserve the
    % variance of _y_, the factor |(1-k)| is neccesary.
    for i = 2:N
       x(i) = a * x(i-1) + rX(i);
       y(i) = a * ((1-k) * y(i-1) + k * x(i-1)) + rY(i) ;
    end
   
    % Correlation coefficient
    C(:,cnt) = correlation(x,y,maxLag);
    cnt = cnt + 1;

end
close(h) % close waitbar

%% Plot results
subplot(3,1,1)
k = linspace(0,1,L); i = 1;
plot(-maxLag:maxLag, C(:,i))
ylim([-1 1]) % limit y-axis to range -1...+1
grid on % show grid lines
xlabel('Lag'), ylabel('Correlation C')
title(sprintf('Coupling k = %0.3f',k(i)))

subplot(3,1,2)
k = linspace(0,1,L); i = 10;
plot(-maxLag:maxLag, C(:,i))
ylim([-1 1]) % limit y-axis to range -1...+1
grid on % show grid lines
xlabel('Lag'), ylabel('Correlation C')
title(sprintf('Coupling k = %0.3f',k(i)))

subplot(3,1,3)
k = linspace(0,1,L); i = 100;
plot(-maxLag:maxLag, C(:,i))
ylim([-1 1]) % limit y-axis to range -1...+1
grid on % show grid lines
xlabel('Lag'), ylabel('Correlation C')
title(sprintf('Coupling k = %0.3f',k(i)))

