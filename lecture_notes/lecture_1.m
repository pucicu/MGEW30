%% Recap Matlab

%% Variables, arrays, cells, structures
x = 1
x = 3;
y = [2 4 8]
y = [2; 4; 8]
z = [1 2 4; 8 16 32; 64 128 256]
z(2,3)
z(:,1)
size(z)

a = {'test', 1}
b{1} = 1:10;
b{2} = 'test'
b{3} = 2

c.flag = 1
c.total = 10:10:100
c.description = 'test'


%% notation (use lower case)
iPosition = 1; % iterator, integer (i, j, k)
x = 0.489; % state, floating point
weight = 0.2; % name should document the meaning
GRAVITY_CONSTANT = 9.81; % constants

%% Arithmetic operators
x + 4

x = [1 2 4]
y = [1; 2; 4]
x * y
y * x
x .* y'

%% Conditions
x = 10
if x < 100
   x
else
   x/2
end


t = 1 | [];          % results in [], so...
if (t) 1, end        % in if ([]), this is false.

if (1 | []) 1, end   % short circuits so condition is true.          

t = [1, 1] | [1, 2, 3];          % error
if ([1, 1] | [1, 2, 3]) 1, end   % OK          


%% Loops
for i = 1:12
  fprintf('%2i: %4i\n',i,2^i)
end

%% Allocating memory
x = zeros(20,1);
y = ones(10,10);

%% Random values
rng(1)
% or randn('seed',1)
x = randn(1,1)
y = rand(10,1)

%% Scripts and functions
% script: write file with extension .m containing
for i = 1:12
  fprintf('%2i: %4i\n',i,2^i)
end

%% function: write file with name power2.m containing
function y = power2(x)
y = x.^2

%% then execute
power2(1:12)


%% Commenting code
% use double %% to structure scripts, use Matlab editor features and
% assign headlines
% use single % for standard comments


%% Plotting
plot(1:10,randn(10,1))
plotyy(1:10,randn(10,1),1:10,100*rand(10,1))

r = randn(10,10);
pcolor(r)

mesh(r)
surf(r)
view(2)

colormap(gray)


%% Layout
% use indentation
% comment from the beginning
% add help text for functions: 
%    1st help line
%    description of function
%    requirements (arguments)
%    example


%% Autoregressive models
% George Udny Yule, 1927, tried to explain and predict the sun spot cycles, second order difference equation
% Gilbert Walker, ~1930, generalised this approach to autoregressive (AR) model for explaining atmospheric circulation: discovery of large-scale circulations, now called Walker circulation

% simple AR model
% x(i) = a * x(i-1) + randn

% in Matlab
clear
N = 1000;
x = zeros(N,1);
r = randn(N,1);
x(1) = r(1);
a = 0.98;
for i = 2:N
   x(i) = a * x(i-1) + r(i);
end

plot(x)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overview Octave
% most functions are the same as in Matlab

%% differences:
t = 1 | [];          % results in [], so...
if (t) 1, end        % in if ([]), this is false.

if (1 | []) 1, end   % short circuits so condition is true.          

t = [1, 1] | [1, 2, 3];          % error
if ([1, 1] | [1, 2, 3]) 1, end   % OK          


%% Random values
randn('seed',1)
x = randn(1,1)
y = rand(10,1)


%% Commenting code
% both % and # can be used


%% Plotting
plot(1:10,randn(10,1))
plotyy(1:10,randn(10,1),1:10,100*rand(10,1))

r = randn(10,10);
pcolor(r)


%% Differences to Matlab
% nested functions
% function handles
% objects, classes 
% no GUIs (but Tcl/Tk and further solutions), no Java
% some differences in graphics
% Matlab file format (.mat) - depends on Octave version 
% commenting character
% short-circuit logical operators
% MEX files not supported
% toolboxes in Matlab and Octave probably not compatible
