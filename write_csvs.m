function otpt=collatergg5(inpt); % $2

listmat=strsplit(ls(strcat(inpt,'/','*mutmat')));
listcov=strsplit(ls(strcat(inpt,'/','*cov')));
listgg=strsplit(ls(strcat(inpt,'/','*mutmatgg')));

% samfile=inpt;
% samidx=findstr(samfile,'/sam/');
% inpt1=samfile(1:samidx+4);


% otpt.sampleid.normat=zeros(21,178); 
% otpt.sampleid.matcum=zeros(21,178); 
% otpt.sampleid.covcum=zeros(1,178);
% otpt.sampleid.ggcum=zeros(1,89);
% otpt.sampleid.wtcum=zeros(1,89);
% otpt.sampleid.normgg=zeros(1,89);

mattmp=zeros(21,178);
covtmp=zeros(1,178);
ggtmp=zeros(1,89);
wttmp=zeros(1,89);

for i=1:length(listmat)-1;
    sampleididx=strfind(listmat{i},'/');
    tmp=listmat{i};
    sampleid=strrep(tmp,'.','_');
    sampleid=strrep(sampleid,'-','_')
    clearvars tmp;

% mattmp=dlmread(listmat{i});
% eval(strcat('otpt.',sampleid,'.matcum=mattmp(1:21,1:178);'));

covtmp1=tdfread(listcov{i});
covtmp=(bin3(covtmp1.cov(:,1)))';
% eval(strcat('otpt.',sampleid,'.covcum=covtmp;'));

ggtmp=dlmread(listgg{i});
% eval(strcat('otpt.',sampleid,'.ggcum=transpose(ggtmp);'));

for k=1:89;
wttmp(1,k)=(covtmp(1,2*k)+covtmp(1,2*k-1))/2;
end
% eval(strcat('otpt.',sampleid,'.wtcum=wttmp;'));

otpt.ggcumall(:,i)=ggtmp';
otpt.wtcumall(:,i)=wttmp';
otpt.labelsall{:,i}=sampleid;

mattmp=zeros(21,178);
covtmp=zeros(1,178);
ggtmp=zeros(1,89);
wttmp=zeros(1,89);
end

for i=1:length(otpt.ggcumall(1,:));
for j=1:89;
otpt.normggallcell{1,i}=otpt.labelsall{i};	
otpt.normggallcell(j+1,i)={otpt.ggcumall(j,i)/otpt.wtcumall(j,i)};
otpt.ggcumallcell{1,i}=otpt.labelsall{i};	
otpt.ggcumallcell(j+1,i)={otpt.ggcumall(j,i)};
otpt.wtcumallcell{1,i}=otpt.labelsall{i};	
otpt.wtcumallcell(j+1,i)={otpt.wtcumall(j,i)};
end
end

% clearvars otpt.ggcumall otpt.labelsall otpt.wtcumall 
% for i=1:length(otpt.ggcumall(1,:));
% subplot(2,5,i);
% scatter(log10(otpt.ggcumall(:,1)),log10(otpt.ggcumall(:,i)));
% title(strrep(otpt.labelsall{1,i},'_',''));
% end
% for i=1:length(otpt.ggcumall(1,:));
% subplot(2,5,i);
% scatter(log10(otpt.normggall(:,1)),log10(otpt.normggall(:,i)));
% title(strrep(otpt.labelsall{1,i},'_',''));
% end

[~,dh_name,~]=fileparts(fileparts(inpt));

 cell2csv(strcat(inpt,dh_name,'_normalised_.csv'),otpt.normggallcell,',');
 cell2csv(strcat(inpt,dh_name,'_counts_.csv'),otpt.ggcumallcell,',');
 cell2csv(strcat(inpt,dh_name,'_coverage_.csv'),otpt.wtcumallcell,',');
disp(strcat('output file : ',inpt,'.csv'))
end
