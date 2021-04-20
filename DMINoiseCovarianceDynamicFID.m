function outputdata=DMINoiseCovarianceDynamicFID(rawdata)
% Read FID data and generate noise covariance matrix
currentfolder=cd;
outputdata.datapath=strcat(currentfolder,'\',rawdata);
[outputdata.data, outputdata.list]=csi_loadData(outputdata.datapath);

%% Parameters
outputdata.Parameters.gyromagneticratio=6.53569*10^6;
outputdata.Parameters.Tesla=7;
outputdata.Parameters.Freq=outputdata.Parameters.gyromagneticratio*outputdata.Parameters.Tesla;
outputdata.Parameters.BW=5000;
outputdata.Parameters.ppmwindow=((outputdata.Parameters.Freq/outputdata.Parameters.BW)^-1)*10^6;
outputdata.Parameters.NP=outputdata.list.F_resolution;
outputdata.Parameters.time= linspace(0,outputdata.Parameters.NP/outputdata.Parameters.BW,outputdata.Parameters.NP);
outputdata.Parameters.zerofill= outputdata.Parameters.NP*2;
outputdata.Parameters.dyn=outputdata.list.number_of_extra_attribute_1_values;

outputdata.xaxis=linspace(-outputdata.Parameters.ppmwindow/2,outputdata.Parameters.ppmwindow/2,outputdata.Parameters.NP);

%% Check signal in tail
outputdata.tails=outputdata.data.raw(0.5*outputdata.Parameters.NP+1:end,:,:,:);
% outputdata.spectra=fftshift(fft(outputdata.tails,[],1),1);
% outputdata.spectratail=fftshift(fft(outputdata.avgtail,[],1),1);
% Tail of the averaged FID is checked if there is still signal found
% No signal found on the last 256 points

%% Concatenate tails of FID to generate a long noise scan
% Aim is to achive at least 20k points
% Each FID contributes 256 points. 80 FID's should be enough
for m=1:size(outputdata.data.raw,2)
outputdata.concatenatedtails(:,m)=reshape(outputdata.tails(:,m,:,:),[numel(outputdata.tails(:,m,:,:)),1]);
end

outputdata.noisecovariance=cov(outputdata.concatenatedtails);
end
