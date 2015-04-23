%% Autoregressive models
clear
randn('seed',1);
N = 1000;
x = zeros(N,1);
y = zeros(N,1);
rX = randn(N,1);
rY = randn(N,1);
x(1) = rX(1);
y(1) = rY(1);
a = 0.98;
k = .60;
for i = 2:N
   x(i) = a * x(i-1) + rX(i);
   y(i) = a * ((1-k) * y(i-1) + k * x(i-1)) + rY(i) ;
end

subplot(3,1,1)
plot(x)
subplot(3,1,2)
plot(y)

subplot(3,1,3)
plot(x,y,'.')
