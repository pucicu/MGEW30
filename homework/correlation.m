function C = correlation(varargin)
% CORRELATION Calculate correlation coefficient between two time series.
%
%   C = CORRELATION(X,Y) calculates the Pearson's correlation coefficient C
%       between the two vectors X and Y.
%
%   C = CORRELATION(A) calculates the Pearson's correlation coefficient C
%       between the first two rows or first two columns of a 2xN or Nx2 
%       matrix A.
%
%   C = CORRELATION(X,Y,TYPE) calculates the correlation coefficient C
%       by using TYPE:
%   Correlation types:
%       'Pearson'   - Pearson's correlaion (default)
%       'Spearman'  - Spearman's rank correlation
%       'Kendall'   - Kendall's tau
%
% Homework 2

% (coding: Norbert Marwan, 5/2015)

%% Check input arguments
% Show error when more input arguments as 3 or no input is given. 
narginchk(1,3)


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
if numel(iDouble) == 1 % matrix given (expected)
    A = double(varargin{iDouble}); % get the numeric input
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
elseif numel(iDouble) == 2 % two time series given as separate vectors
    x = double(varargin{iDouble(1)});
    y = double(varargin{iDouble(2)});
    % error if X and Y are of different length
    if length(x) ~= length(y)
        error('Time series X and Y must be of equal length.')
    end
else
    % this should not happen
    error('No valid numerical input given.')
end


%% Calculate correlation
switch(type)
  case 'pearson'
    C = pearson(x,y);
    
  case 'spearman'
    warning('Spearman''s rank correlation not yet implemented.')
    C = NaN;
    
  case 'kendall'
    warning('Kendall''s tau not yet implemented.')
    C = NaN;
    
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

