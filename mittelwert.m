function m = mittelwert(varargin)
% MITTELWERT   Berechnet den Mittelwert von Vektoren und Matrizen
% M = MITTELWERT(X,Y) berechnet den Mittlwert M aus dem Vektor X
% oder Matrix X. X kann Spalten- oder Zeilenvektor sein, oder sogar eine
% Matrix. Das zweite Argument ist zur Demonstration der Matlab-Variablen
% VARARGIN. Y wird einfach mit dem Mittelwert multipliziert, wobei die
% Angabe von Y optional ist; wird es nicht angegeben, ist Y = 1.

% N. Marwan, 23. April 2015
%% Eingabe-Check
% testet wieviele Eingabeargumente angegeben wurden.
narginchk(1,2)

%% Eingabewerte einlesen
% varargin ermöglicht die flexible Übergabe von Eingabewerten
% das erste Argument ist die Zeitreihe aus der der Mittelwert
% berechnet wird
x = varargin{1};

%%
% falls noch ein zweites Argument angegeben wurde, soll dieses in die
% Variable y übergeben werden, ansonsten ist y = 1
if nargin == 2
    y = varargin{2};
else
    y = 1;
end

%% Dimension x
% N = Anzahl Zeilen, M = Anzahl Spalten
[N M] = size(x);

%% Initialisiere Ausgabevektor
% je nachdem, ob Vektor oder Matrix, wird die Ausgabe ein Skalar 
% oder ein Vektor sein
m = zeros(1,M);

%% Berechnung
% Aufsummieren der Werte in jeder Spalte separat
for j = 1:M
    for i = 1:N
       m(j) = m(j) + x(i,j);
    end
end

% normiere durch Länge
m = m/N;

% skaliere Mittelwert mit y
m = y * m;
