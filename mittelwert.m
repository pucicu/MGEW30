function m = mittelwert(varargin)
% MITTELWERT ?

narginchk(1,2)

x = varargin{1};

if nargin == 2
    y = varargin{2};
else
    y = 1;
end
[N M] = size(x);
m = zeros(1,M);
for j = 1:M
    for i = 1:N
       m(j) = m(j) + x(i,j);
    end
end

m = m/N;

m = y * m;
