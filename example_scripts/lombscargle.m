function [Pxx Fxx Hxx] = lombscargle(t,x)
% LOMBSCARGLE   Computes the Lomb-Scargle periodogram.
%    [PXX FXX HXX] = LOMBSCARGLE(T,X) computes the Lomb-Scargle periodogram
%    of unevenly spaced data in vector X. The vector T contains the
%    unevenly spaced time scale. The results are the periodogram PXX,
%    the corresponding frequencies are in FXX, and the 
%    corresponding 95% false alarm levels HXX.
%
%    Reference:
%    Trauth, M. H.: MATLAB Recipes for Earth Sciences, Springer 2007.


t = t(:);
y = x(:);

int = mean(diff(t));
ofac = 4; hifac = 1;
f = ((2*int)^(-1))/(length(y)*ofac): ...
    ((2*int)^(-1))/(length(y)*ofac): ...
    hifac*(2*int)^(-1);

y = y - mean(y);
vy = var(y);

for k = 1:length(f)
    wrun = 2*pi*f(k);
    fact1 = wrun*t - atan2(sum(sin(2*wrun*t)),sum(cos(2*wrun*t)))/2;
    
    px(k) = 1/(2*vy) * ...
       ((sum(y.*cos(fact1))).^2) /(sum((cos(fact1)).^2)) + ...
       ((sum(y.*sin(fact1))).^2) /(sum((sin(fact1)).^2));
end

prob = 1-(1-exp(-px)).^length(y);

m = floor(0.5*ofac*hifac*length(y));
effm = 2*m/ofac;
signif = 0.95;
levels = log((1-signif.^(1/effm)).^(-1));

plot(f,px)
hold on
for k = 1:length(signif)
    line(f,levels(:,k)*ones(size(f)),'LineStyle','--')
end
xlabel('Frequency')
ylabel('Power')
title('Lomb-Scargle Powerspectrum')
hold off

if nargout > 0
    Pxx = px;
end
if nargout > 1
    Fxx = f;
end
if nargout > 2
    Hxx = levels;
end
    
