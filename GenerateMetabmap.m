function MetabMap=GenerateMetabmap(AMARESoutput,metabnum)
Metabnames={'Lipid';
%     'Glx';
    'Glucose';
    'Water'};
disp(strcat('Generating metabolite map: ',cellstr(Metabnames{metabnum})))
a=cell2mat(AMARESoutput);
b=zeros(numel(a),3);
for m=1:numel(a)
    b(m,:)=a(m).amplitude;
end
c=b(:,metabnum);
MetabMap=reshape(c,size(AMARESoutput));
