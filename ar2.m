%% Coupled Autoregressive System
% Example of two coupled systems as a test model for homework.
% It is also used to show the _publish_ functionality of MATLAB.

%% Initialisation
% Some first settings.
clear
randn('seed',1); % reset the random number generator to ensure reproducable random numbers
N = 1000; % length of the time series
x = zeros(N,1); % initialize vector with N times the value of zero
y = zeros(N,1); % initialize vector with N times the value of zero
rX = randn(N,1); % create random numbers from Gaussian distribution
rY = randn(N,1); % create random numbers from Gaussian distribution
x(1) = rX(1); % initialze the first value in the time series with a random number
y(1) = rY(1); % initialze the first value in the time series with a random number


%% Parameter
% These parameters should be changed in the home work in order to
% study the correlation vs. a change of _a_ resp. _k_.
a = 0.98; % AR coefficient (for both systems the same)
k = .10; % coupling constant

%% Calculate the coupled system
% The system _x_ is a standard AR(1), the system _y_ is coupled with
% system _x_ by a coupling constant _k_. In order to preserve the
% variance of _y_, the factor |(1-k)| is neccesary.
for i = 2:N
   x(i) = a * x(i-1) + rX(i);
   y(i) = a * ((1-k) * y(i-1) + k * x(i-1)) + rY(i) ;
end

%% Plot results
% Use |subplot| to have separate panels in the figure window.
subplot(3,1,1)
plot(x)
ylabel('System x')

subplot(3,1,2)
plot(y)
ylabel('System y')

subplot(3,1,3)
plot(x,y,'.')
xlabel('System x'), ylabel('System y')
