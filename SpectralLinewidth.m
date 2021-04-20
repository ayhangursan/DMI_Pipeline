function [T2starmap, Hzmap, ppmmap]=SpectralLinewidth(spectradata,xaxis,Parameters)

% poolobj = gcp('nocreate'); % Start parallel pools

%% Shorten vector size
waterwindow=find(xaxis>4 & xaxis<5.4);

xaxisshortened=xaxis(waterwindow);
spectradatashortened=real(spectradata(waterwindow,:));


%%
Hzmap=zeros(Parameters.CSIdims);
ppmmap=zeros(Parameters.CSIdims);
interp_factor = 10^-2;
xaxisint= interp1(1:size(xaxisshortened,2),xaxisshortened,1:interp_factor:size(xaxisshortened,2));
%%
for m=1:prod(Parameters.CSIdims)
    spectraint=real(interp1(xaxisshortened,spectradatashortened(:,m),xaxisint));
    [max_value, max_ind] = max((spectraint));
     mean_value = mean((spectraint));
if max_value > mean_value
    max_half = max_value/2;
    [minValueR,closestIndexR] = min(abs(spectraint(max_ind:end)-max_half)); % Right side of the fit
    closestIndexR=closestIndexR+max_ind;
    [minValueL,closestIndexL] = min(abs(spectraint(1:max_ind)-max_half)); % Left side of the fit
    if closestIndexR>numel(xaxisint)
       ppmmap(m)=0;
    else
    ppmmap(m)=abs(xaxisint(closestIndexL)-xaxisint(closestIndexR));
    end
end
    clear spectraint;
end
%%

Hzmap=ppmmap*Parameters.Freq/10^6;
T2starmap=1000./(pi*Hzmap);% in ms
% delete(poolobj); %Shut parallel pools
% T2starmap=Hzmap

end