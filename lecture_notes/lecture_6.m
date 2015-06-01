%% Functions

% inline functions and function handles
f = inline('1-a./x.^2','a','x')
f = @(a,x)[1-a./x.^2]

x = linspace(0,10,100);
plot(x,f(2,x))

%% Root finding

%% bisection method
f = inline('1-1./x.^2','x')

a = 0;
b = 20;
for i = 1:10
    x(i) = (a + b)/2;
    if f(x(i)) > 0
        b = x(i);
    elseif f(x(i)) < 0
        a = x(i);
    end
    x(i)
end
plot(x)


%% Newton root finding
f = inline('1-1./x.^2','x')
f1 = inline('2./x.^3','x');
x = 1.5;
for i = 1:10
    x(i+1) = x(i) - f(x(i)) / f1(x(i));
end

plot(x)


%% Matlab version
f = inline('1-2./x.^2','x')
fzero(f,2)



%% polynoms
p = [1 0 -12 0 .5]
%^p = [1 -6 11 -6];
x = linspace(-5,5,100);
y = polyval(p,x);
plot(x,y)

% roots
roots(p)

% integral and derivative of polynoms
p = [3 1 0 0];
polyint(p)
polyder(p)

% empirical derivation from time series
diff(x)

% fitting a polynom to data
x = load('example.dat');
plot(x(:,1),x(:,2))
xlabel('Depth')
ylabel('Age')

p = polyfit(x(:,1),x(:,2),3);

plot(x(:,1),x(:,2),x(:,1),polyval(p,x(:,1)))
xlabel('Depth')
ylabel('Age')


%% iterative processes/ fixed point iterations

f = inline('cos(x)','x')

f = inline('1-cos(x)','x')


%% sedimentation example

h = 0;
r = 5; % sedimentation rate, mm/yr
N = 10; % 10 years

for i = 2:N
    h(i) = h(i-1) + r + 2*rand;
end
plot(h)

% time-dependent rate
r = 5*cos(linspace(0,pi/4,N));
for i = 2:N
    h(i) = h(i-1) + r(i) + 5*rand;
end
plot(h)

% plot with some random colours
C = jet(N);
C = C(randperm(N),:);

clf
hold on
for i = 1:N-1
    fill([0 1 1 0 0],[h(i) h(i) h(i+1) h(i+1) h(i)],C(i,:))
end







%% population dynamics

%% Malthusian growth model
% q = 1-d+b (d = death rate, b = birth rate)

x0 = 5;
x = x0;
q = 1.1;
for i = 1:500
   x(i+1) = q * x(i);
end

% or for an initial state x0
x = inline('q.^n * x0','q','x0','n');
n = linspace(0,5,50);
plot(n,x(x0,q,n))

% food production
k = 500;
y = inline('k.*n','k','n');
plot(n,x(x0,q,n),n,y(k,n))

% leads to Malthusian catastrophe!


%% crossing point with Newton scheme
f = inline('q.^n * x0 - k * n','q','x0','k','n');
f1 = inline('n.*q.^(n-1) * x0 - k','q','x0','k','n');
x = x0;
for i = 1:1000
    x(i+1) = x(i) - f(x0,q,k,x(i)) / f1(x0,q,k,x(i));
end

plot(x)


%% growth model with feedback (logistic map)
% q is not constant but depends on x, e.g., q = (1-x)*r

x(i+1) = r * x(i) * (1 - x(i))

f = inline('r * x .* (1 - x)','r','x')

f(3,.2)
x = .2;
for i = 1:1000
    x(i+1) = f(3.8,x(i));
end
plot(x)


% now the dynamics depends on r (test for different r)


%% Integration

%% midpoint integration

% limits a, b
% M number of subintervals
% H length of subinterval

% simple example:
f = inline('-x+1','x')
a=0;
b=1;
% other example:
f = inline('cos(x)','x')
a = -pi/2;
b = pi/2;

% integration
M = 10;
H = (b-a)/M;
x = linspace(a+H/2,b-H/2,M); 
y = f(x);
Imp = H*sum(y)

%% trapezoidal integration

f = inline('cos(x)','x')
a = -pi/2;
b = pi/2;
M = 1000;
H = (b-a)/M;
x1 = linspace(a,b-H,M); 
x2 = linspace(a+H,b,M); 
y = (f(x1) + f(x2));
Imp = sum(y) * H/2

%% Simpson integration
f = inline('cos(x)','x')
a = -pi/2;
b = pi/2;
M = 1;
H = (b-a)/M;
x1 = linspace(a,b-H,M); 
x3 = linspace(a+H,b,M); 
x2 = (x1+x3)/2; 
y = (f(x1) + 4*f(x2) + f(x3));
Imp = sum(y) * H/6



%% Matlab version

% (works only with function handles)
f = @(x)[cos(x)]
a = -pi/2;
b = pi/2;
integral(f,a,b)

%% Integration on time series
trapz
cumtrapz

a = -pi/2;
b = pi/2;
x = linspace(a,b,100);
y = f(x);

trapz(y) * mean(diff(x))
