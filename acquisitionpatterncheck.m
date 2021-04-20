function patternarray=acquisitionpatterncheck(rawdata)
dims=size(rawdata);
weightedpattern=find(squeeze(sum(abs(rawdata),1)));
B=zeros(dims(2:end));
B(weightedpattern)=1;
patternarray=sum(B,ndims(B));
