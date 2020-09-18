if ~exist('pathName','var')
    pathName=pwd;
    if ~strcmpi(pathName(end),filesep)
        pathName=[pathName,filesep];
    end
end
[fileName, pathName]=uigetfile([pathName,'*.csv'],'Get main data file.');
csv490File=[pathName,fileName];
[fileName, pathName]=uigetfile([pathName,'*.csv'],'Get control data file.');
csv405File=[pathName,fileName];

if ~isempty(regexp(fileName,'chl(\d)'));
    findString='chl(\d)';
    replaceString='avg';
else
    findString='.csv';
    replaceString='_avg.csv';
end
outFileName=regexprep(fileName,findString,replaceString);
outFileName=[pathName,outFileName];


mat490=csvread(csv490File,1,3);
mat405=csvread(csv405File,1,3);

mat490=mat490(:,find(any(mat490,1)));
mat405=mat405(:,find(any(mat405,1)));

if ~all(size(mat490)==size(mat405))
    error('Files must have the same number of entries.');
end


timeVals=mat490(:,1);

mat490=mat490(:,5:end);
mat405=mat405(:,5:end);

mat490Avg=mean(mat490,2);
mat405Avg=mean(mat405,2);

outMat=[timeVals,mat490Avg,mat405Avg];

csvwrite(outFileName,outMat);