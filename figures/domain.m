datadir = 'C:\Users\labpatron\Desktop\CODE\CODE\no_evap_3\';

filePattern = fullfile(datadir, 'steak_default.0*.mat');
matFiles = dir(filePattern);
names = matFiles.name;

ddomain = NaN(length(matFiles),2,2);
figure(1);
hold on;

for k = 1:length(matFiles)
        load(sprintf('%s%s',datadir,matFiles(k).name));
     ddomain(k,1,1)= min(S.T(2,2:end-1));
     ddomain(k,1,2)= max(S.T(2,2:end-1));
     ddomain(k,2,1)= min(S.phi(2,2:end-1));
     ddomain(k,2,2)= max(S.phi(2,2:end-1));
     rectangle('Position',[ddomain(k,1,1) ddomain(k,2,1) ddomain(k,1,2)-ddomain(k,1,1) ddomain(k,2,2)-ddomain(k,2,1)])
end