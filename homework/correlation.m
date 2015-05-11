function varargout = correlation(varargin)
% CORRELATION Calculate correlation coefficient between two time series.
%
%   C = CORRELATION(X,Y) calculates the Pearson's correlation coefficient C
%       between the two vectors X and Y.
%
%   C = CORRELATION(A) calculates the Pearson's correlation coefficient C
%       between the first two rows or first two columns of a 2xN or Nx2 
%       matrix A.
%
%   C = CORRELATION(A,MAXLAG...) calculates the cross correlation coefficient C
%       over a range of lags up to a maximal lag given by the scalar MAXLAG.
%
%   [C L] = CORRELATION(A,MAXLAG,...) returns the vector of lags L.
%
%   C = CORRELATION(X,Y,...,TYPE) calculates the correlation coefficient C
%       by using TYPE:
%   Correlation types:
%       'Pearson'   - Pearson's correlaion (default)
%       'Spearman'  - Spearman's rank correlation
%       'Kendall'   - Kendall's tau
%
%   Example:  Calculate Spearman correlation between two random vectors.
%       x = rand(100,1);                    % time series X
%       y = rand(100,1);                    % time series Y
%       C = correlation(x,y,'spearman');    % compute Spearman correlation
%
%       A = randn(1000,2);                  % two time series in A
%       [C L] = correlation(A,5,'kendall'); % compute Kendall's tau cross correlation
%

% Homework 3
% (coding: Norbert Marwan, 5/2015)

%% Check input/output arguments
% Show error when more input arguments than 3 or no input is given. 
narginchk(1,4)

% Show error when more output arguments than 2. 
nargoutchk(0,2)


%% Read input variables
% Find index of number inputs (the data, either X, Y, or As).
iDouble=find(cellfun('isclass',varargin,'double'));
% Find index of string inputs (the TYPE).
iChar=find(cellfun('isclass',varargin,'char'));

% Check for variable TYPE
type = 'pearson';
if ~isempty(iChar)
    type = lower(varargin{iChar}); % get char input and convert to lowercase
    if ~ismember(type,{'pearson','spearman','kendall'}) 
    % check whether the input corresponds to one of the three allowed TYPEs.
    % If not set to default and show a warning message.
        type = 'pearson';
        warning('Could not recognize correlation type. Set to default (Pearson).')
    end
end

% Check for the numerical input data
maxlag = 0;
if numel(iDouble) == 1 | (numel(iDouble) == 2 & numel(double(varargin{iDouble(2)})) == 1) % matrix given (expected)
    A = double(varargin{iDouble(1)}); % get the numeric input
    N = size(A); % size of the input
    if min(size(A)) ~= 2 % matrix must be 2-columns or 2-rows, else error
        error('At least one dimension of A must be 2.')
    end
    % now read data from matrix, depending whether 2-columns or 2-rows
    if N(2) == 2 % 2-columns
        x = A(:,1);
        y = A(:,2);
    else % 2-rows, transpose to ensure column vectors
        x = A(1,:)';
        y = A(2,:)';
    end
    % if MAXLAG is provided, it should be the 2nd element
    if numel(iDouble) == 2 & numel(double(varargin{iDouble(2)})) == 1
        maxlag = double(varargin{iDouble(2)});
    end
    
elseif numel(iDouble) >= 2 & numel(double(varargin{iDouble(2)})) ~= 1 % two time series given as separate vectors
    x = double(varargin{iDouble(1)});
    y = double(varargin{iDouble(2)});
    % error if X and Y are of different length
    if length(x) ~= length(y)
        error('Time series X and Y must be of equal length.')
    end
    if numel(iDouble) == 3 % the 3rd numeric argument is the MAXLAG
        maxlag = double(varargin{iDouble(3)});
    end
    
else
    % this should not happen
    error('No valid numerical input given.')
end

% check whether the maxlag is negative
if maxlag < 0
    warning('MAXLAG cannot be negative.')
    maxlag = abs(maxlag);
end


%% Calculate correlation
C = zeros(2*maxlag + 1,1);

% correlation for shifted time series
for i = 1:maxlag
    
    % shift y in negative direction
    xCutted = x((1+i):end);
    yShifted = y(1:end-i);
    
    switch(type)
      case 'pearson'
        C(maxlag-i+1) = pearson(xCutted,yShifted);

      case 'spearman'
        C(maxlag-i+1) = spearman(xCutted,yShifted);

      case 'kendall'
        C(maxlag-i+1) = kendall(xCutted,yShifted);

    end

    % shift y in positive direction
    xCutted = x(1:end-i);
    yShifted = y((1+i):end);
    switch(type)
      case 'pearson'
        C(i+maxlag+1) = pearson(xCutted,yShifted);

      case 'spearman'
        C(i+maxlag+1) = spearman(xCutted,yShifted);

      case 'kendall'
        C(i+maxlag+1) = kendall(xCutted,yShifted);
    end

end

% finally the correlation without lag
switch(type)
  case 'pearson'
    C(maxlag+1) = pearson(x,y);

  case 'spearman'
    C(maxlag+1) = spearman(x,y);

  case 'kendall'
    C(maxlag+1) = kendall(x,y);

end

%% put results into output variables
varargout{1} = C;
% if two output variables are available
if nargout == 2
    varargout{2} = (-maxlag:maxlag)';
end


%% Section with some helper functions
function C = pearson(x,y)
% Pearson's correlation coefficient

    L = length(x); % length of the time series
    % mean
    xMean = sum(x)/L; % mean for X
    yMean = sum(y)/L; % mean for Y
    % stddev
    xStd = sqrt(sum((x - xMean).^2) / (L-1)); % standard deviation for X
    yStd = sqrt(sum((y - yMean).^2) / (L-1)); % standard deviation for Y
    % covariance
    cv = sum((x - xMean) .* (y - yMean))/(L-1);
    % correlation
    C = cv / (xStd * yStd); % Pearson correlation is the normalised covariance


function C = spearman(x,y)
% Spearman's correlation coefficient, only for time series with unique values

    L = length(x); % length of the time series
    % assign ranks
    [~,~,xrank] = unique(x);
    [~,~,yrank] = unique(y);
    
    % check for consistency
    if length(xrank) ~= length(yrank)
        error('Time series do not contain unique values.')
    end    
    
    % check for unique values
    if L ~= numel(unique(x)) | L ~= numel(unique(x))
        error('Time series do not contain unique values.')
    end
    
    % correlation
    sumdsq = sum((xrank-yrank).^2);
    C = 1-6*sumdsq/(L*(L^2-1)); %this is the Rho of the spearman correlation


function C = kendall(x,y)
% Kendall's correlation coefficient, only for time series with unique values

    L = length(x); %length of the time series
    
    % check for unique values
    if L ~= numel(unique(x)) | L ~= numel(unique(x))
        error('Time series do not contain unique values.')
    end
    
    
    % sort first time series and reorder the second in the same order
    [xsort,idx] = sort(x);
    ysort = y(idx);
    
    numcon = 0; % concordant pairs

    % count concordant pairs, because x is sorted, we just need to check y
    for i = 2:L
        %numcon = numcon + sum(0.5* (1+sign(ysort(i:end)-ysort(i-1))));
        numcon = numcon + sum(ysort(i:end)>ysort(i-1));
    end

    numdis=(L^2-L)/2-numcon;
    
    % correlation
    C = 2*(numcon-numdis)/(L^2-L);


