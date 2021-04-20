function DMIResults=SplittingDetectionAmares(inputfids,Parameters,xaxis,waterloc)
defaultfolder=cd;
AMARESfolder='C:\Users\agursan.DS\Desktop\DMI\Matlab_codes\OXSA-master\main\+AMARES';
DMIpar.samples=Parameters.NP;
DMIpar.imagingFrequency=Parameters.Freq/(10^6);
DMIpar.timeAxis=Parameters.time;
DMIpar.dwellTime=1/(Parameters.BW);
DMIpar.ppmAxis=xaxis-4.7;
DMIpar.offset=-4.7; % This may be modified in the future to include frequency shift based on B0 imperfections
DMIpar.beginTime=Parameters.TE;disp(strcat('TE=',num2str(Parameters.TE*1000),'ms for AMARES'));
% DMIpar.beginTime=0.00195;
matrixsize=size(inputfids);
CSIdims=matrixsize(2:end);
% Load prior knowledge file
cd('C:\Users\agursan.DS\surfdrive\Deuterium_CSI\Matlab_codes\AMARESforDMI\'); % change folder to load prior knowledge
pk = PK_2H_7T_QuadrupolarSplitting;disp('Lorantzian lineshapes in Brain.') % Which prior knowledge will be used is to be decided

pkglobal=pk;

cd(AMARESfolder)

% for i=1:length(pk.initialValues)
%     pk.initialValues(i).chemShift = (pk.initialValues(i).chemShift + DMIpar.offset);
%     pk.bounds(i).chemShift = (pk.bounds(i).chemShift + DMIpar.offset);
% end

for fidnum=1:numel(inputfids)/DMIpar.samples
    for i=1:length(pk.initialValues)
        
        DMIpar.offset=-(4.7+(4.7-waterloc(fidnum))); % This may be modified in the future to include frequency shift based on B0 imperfections
        pk.initialValues(i).chemShift = (pkglobal.initialValues(i).chemShift + DMIpar.offset);
                pk.bounds(i).chemShift = (pkglobal.bounds(i).chemShift + DMIpar.offset);
    end
    
    [fitResults, fitStatus, ~, CRBResults] = AMARES.amaresFit(inputfids(:,fidnum), DMIpar, pk, 0);
    [row,col,slice]=ind2sub(CSIdims,fidnum);
    
    DMIResults{row,col,slice}=fitResults;
    cd(defaultfolder)
end

disp('AMARES fitting finished.')
end