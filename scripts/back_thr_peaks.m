function thr = back_thr_peaks(x,R,G,B)
%BACK_REM_PEAKS Estimatre threshold on each color channel to remove
%background

thr = zeros(1,3);
thr(1) = find_thr(x,R);
thr(2) = find_thr(x,G);
thr(3) = find_thr(x,B);

end

function thr = find_thr(x,y)

% remove no signal at the beginning or the end
ind = find(y>0);
ind = ind(1):ind(end);
x = x(ind);
y = y(ind);

% smooth signal to reduce compression artifacts
% figure; hold on; plot(x,y);plot(x,smoothdata(y,'gaussian',9))
y = smoothdata(y,'gaussian',9);

% find peaks
[xpeak,~,minima] = peakdetect(x,y,0,1.25,0.001);

% change from indices to x-values
if ~isempty(xpeak); [~,xpeak] = ismember(xpeak,x); end

% select last minimum as threshold
if isempty(minima)     %only 1 peak
    corr = true;
else
    thr = x(minima(end));
    
    ind = y>0;
    if thr < median(x(ind))
        corr = true;
    else
        corr = false;
        
    end
end

%probably only signal peak
if x(xpeak(end)) < 225
    xpeak = [xpeak(end),length(x)];
    ind = xpeak(1):xpeak(2);
    [~,ind2] = min(y(ind));
    thr = x(ind(ind2));
    corr = false;
end


% if needs correction, use double knee method
if corr
    if isempty(minima)
        minima = 1;
    end
    ind = minima(end):xpeak(end);
    
    ind2 = find_inflection(x(ind),y(ind));
    ind = minima:ind(ind2);
    ind2 = find_inflection(x(ind),y(ind));
    thr = x(ind(ind2));
    
    %   double  Otsu
    %   tmp = max(x)*otsuthresh(y);
    %   thr = round(tmp + (max(x)-tmp)*otsuthresh(y(tmp:end)));
end
end

function [x,y,minima] = peakdetect(x,y,thr,cond1,cond2)
% PEAKDETECT(X,Y,THR,COND1,COND2)
%
% This function founds peaks in a mass spectrum spectrum
% using first derivative and then reduces its number using
% small peaks reduction, low-intensity peaks reduction and/
% noise peaks reduction.
%
% IN:
%   x - x data [nx1]
%   y - y data [nx1]
%   thr - intensity threshold value [0.01]
%   cond1 - conditon for small peaks reduction [1.25]
%   cond2 - condition for noise peaks reduction [0.001]
%
% Author:
% Michal Marczyk
% michal_marczyk@op.pl

if nargin < 2;
    disp('Not enough input parameters!');
    x = [];
    y = [];
    minima = [];
    return;
end

if nargin < 5
    disp('Default parameter values will be used.')
    thr = 0.01;
    cond1 = 1.25;
    cond2 = 0.001;
end

x = x(:);
y = y(:);
thr = thr * max(y);    % threshold for low intensities peaks removal
sign = diff(y)>0;      % finds all extrema
maxima = find(diff(sign)<0) + 1;                %all maxima
minima = [find(diff(sign)>0)+1 ;length(y)];     %all minima
if length(maxima) < 1
    [~,maxima] = max(y);
end
if minima(1) > maxima(1)
    minima = [1 ;minima];
end
if maxima(end) == length(y)
    maxima(end) = maxima(end)-1;
end

% small peaks removal by slopes
ys = y;
ys(ys<1) = 1;
count = 1;
peak = zeros(1,length(maxima));
for i=1:length(maxima)
    y_check = y(maxima(i));
    if abs(y_check/ys(minima(i)))>cond1 || abs(y_check/ys(minima(i+1)))>cond1
        peak(count) = maxima(i);
        count = count + 1;
    end
end
peak(count:end) = [];

% low-intensity peaks removal
peak = peak(y(peak) > thr);

i = 1;
%noise peaks removal
while i<length(peak)
    match = check([x(peak(i)) x(peak(i+1))],[y(peak(i)) y(peak(i+1))],cond1,cond2);
    startpos = i;
    while match == 1
        i = i+1;
        if i<length(peak)
            match = check([x(peak(startpos)) x(peak(i+1))],[y(peak(i)) y(peak(i+1))],cond1,cond2);
        else
            break
        end
    end
    endpos = i;
    if endpos-startpos > 0
        sign = diff(y(peak(startpos)-1:peak(endpos)+1))>0;
        maxima = find(diff(sign)<0) + 1;
        [~, ind] = max(y(peak(startpos)+maxima(:)-2));
        peak(endpos) = peak(startpos) - 2 + maxima(ind);
        peak(startpos:endpos-1) = 0;
    end
    i = endpos + 1;
end
peak(peak == 0) = [];

x = x(peak);

% find minima between peaks
n = length(x);
minima = zeros(n-1,1);
for a=1:(n-1)
    [~,ind] = min(y(peak(a):peak(a+1)));
    minima(a) = peak(a) + ind;
end
y = y(peak);

end

function match = check(tempx,tempint,cond1,cond2)
% match test under 2 conditions
match = 0;
if tempint(1)/tempint(2) < cond1 && tempint(1)/tempint(2) >  1/cond1
    if abs(tempx(1) - tempx(2)) < cond2*tempx(1)
        match = 1;
    end
end
end