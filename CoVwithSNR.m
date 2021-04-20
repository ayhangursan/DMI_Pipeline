function []=CoVwithSNR(SNRmap,CoVMap)
% SNRmap=V9158.Scan1.SNR.Roemer;
MatrixDims=size(SNRmap);
% CoVMap=V9158.CoV;

[Max10SNR, Max10SNRLoc]=maxk(SNRmap(:),10);
[x, y, z]=ind2sub(MatrixDims,Max10SNRLoc);

for m=1:10
    CoVonHighestSNR10(m)=CoVMap(x(m),y(m),z(m));
end

[Max50SNR, Max50SNRLoc]=maxk(SNRmap(:),50);
[x, y, z]=ind2sub(MatrixDims,Max50SNRLoc);

for m=1:50
    CoVonHighestSNR50(m)=CoVMap(x(m),y(m),z(m));
end


% SNRmap(find(SNRmap<66.5 & SNRmap>53.4))
% CoVonmatchingSNR10=CoVMap(find(SNRmap<73 & SNRmap>70.5));
% CoVonmatchingSNR50=CoVMap(find(SNRmap<66.5 & SNRmap>53.4));


disp('Coefficient of variation in 10 voxels with highest SNR:');
disp(strcat(num2str(mean(CoVonHighestSNR10)),' std:',num2str(std(CoVonHighestSNR10))))
disp(strcat('SNR values:',num2str(mean(Max10SNR)),' std:',num2str(std(Max10SNR))))

disp('Coefficient of variation in 50 voxels with highest SNR:');
disp(strcat(num2str(mean(CoVonHighestSNR50)),' std:',num2str(std(CoVonHighestSNR50))))
disp(strcat('SNR values:',num2str(mean(Max50SNR)),' std:',num2str(std(Max50SNR))))

% disp('Phantom data. Coefficient of variation in 10 voxels with matching SNR:');
% disp(strcat(num2str(mean(CoVonmatchingSNR10)),' std:',num2str(std(CoVonmatchingSNR10))))
% disp(strcat('SNR values:',num2str(mean(SNRmap(find(SNRmap<73 & SNRmap>70.5)))),' std:',num2str(std(SNRmap(find(SNRmap<73 & SNRmap>70.5))))))
% 
% disp('Phantom data. Coefficient of variation in 50 voxels with matching SNR:');
% disp(strcat(num2str(mean(CoVonmatchingSNR50)),' std:',num2str(std(CoVonmatchingSNR50))))
% disp(strcat('SNR values:',num2str(mean(SNRmap(find(SNRmap<66.5 & SNRmap>53.4)))),' std:',num2str(std(SNRmap(find(SNRmap<66.5 & SNRmap>53.4))))))
end
