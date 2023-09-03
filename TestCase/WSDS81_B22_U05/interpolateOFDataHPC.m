function []=interpolateOFDataHPC(intfOption,sampleNbr,subDomainNbr)


caseInfo.intfOption = intfOption;
caseInfo.sampleNbr = sampleNbr;
caseInfo.subDomainNbr = subDomainNbr;

switch caseInfo.intfOption
    case 1
        caseInfo.intfVar = 'HeavisideInf';
        caseInfo.intfLim = 0.5;
    case 2
        caseInfo.intfVar = 'psi';
        caseInfo.intfLim = 0;
    otherwise
        caseInfo.intfVar = 'alpha.water';
        caseInfo.intfLim = 0.5;
end

caseInfo.CAVar = 'thetaDeg';

caseInfo.fullPathName=pwd;
[caseInfo]=getDropProperties(caseInfo);
[dirInfo]=findDirInfo(caseInfo);
[caseInfo,dirInfo]=getDMax(caseInfo,dirInfo);
[caseInfo,dirInfo]=getDSample(caseInfo,dirInfo,'sigma',1,caseInfo.sampleNbr);

tickzSigma=cell2table(cell(length(caseInfo.nTs)+1,2),'VariableNames',{'time','sigma'});
for k=1:height(tickzSigma)
    if (k==1)
        kMod=caseInfo.nTi;
    else
        kMod=caseInfo.nTs(k-1);
    end
    tickzSigma.time{k}=(dirInfo.time(kMod)-dirInfo.time(caseInfo.nT0))*1000;
    tickzSigma.sigma{k}=cat(2,dirInfo.intf{kMod}/caseInfo.dropRadius,dirInfo.sigma{kMod}*1000);
end

save(['tikz','.mat'],'caseInfo','tickzSigma','-v7.3') %No dirInfo to save opening file time 
%save(['tikz','.mat'],'caseInfo','dirInfo','tickzSigma','-v7.3')

end

%%
function [caseInfo,dirInfo]=getDSample(caseInfo,dirInfo,varName,varDim,nS)

Ds=linspace(0,dirInfo.D{caseInfo.nTmax},nS);
nTs=zeros(size(Ds));
nTs(1)=caseInfo.nTmin;
nTs(nS)=caseInfo.nTmax;

for k=2:nS-1
    
    IMold=0;
    DSwitch=true;
    while DSwitch
        I=find(cellfun(@isempty, table2cell(dirInfo(:,'D')))==0);
        I=I(and((I>=caseInfo.nTmin),(I<=caseInfo.nTmax))); 
        D=cat(1,dirInfo.D{I});
        
        f=find(Ds(k)-D>=0,1,'last');
        if isempty(f)
            IL=caseInfo.nTmin;
            IH=caseInfo.nTmax;
        else
            IL=I(f);
            if (f+1<length(I))
                IH=I(f+1);
            else
                IH=caseInfo.nTmax;
            end
        end
 
        switch (IH-IL)
            case 0
                IM=IL;
                DSwitch=false;
            case 1
                if (Ds(k)-dirInfo.D{IL}<dirInfo.D{IH}-Ds(k))
                    IM=IL;
                else
                    IM=IH;
                end
                DSwitch=false;
            otherwise
                IM=floor(IL+(Ds(k)-dirInfo.D{IL})*(IH-IL)/(dirInfo.D{IH}-dirInfo.D{IL}));
                if (IM==IMold)
                    if isempty(dirInfo.D{IM+1})
                        timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{IM+1}];
                        geoInfo=getGeo(timePathName);
                        [inftC,inftMax]=getIntfC(caseInfo,geoInfo,timePathName);
                        dirInfo.geo{IM+1}=geoInfo;
                        dirInfo.intf{IM+1}=inftC;
                        dirInfo.D{IM+1}=inftMax;
                    end
                    if isempty(dirInfo.D{IM-1})
                        timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{IM-1}];
                        geoInfo=getGeo(timePathName);
                        [inftC,inftMax]=getIntfC(caseInfo,geoInfo,timePathName);
                        dirInfo.geo{IM-1}=geoInfo;
                        dirInfo.intf{IM-1}=inftC;
                        dirInfo.D{IM-1}=inftMax;
                    end
                end
        end
         
        if isempty(dirInfo.D{IM})
            timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{IM}];
            geoInfo=getGeo(timePathName);
            [inftC,inftMax]=getIntfC(caseInfo,geoInfo,timePathName);
            dirInfo.geo{IM}=geoInfo;
            dirInfo.intf{IM}=inftC;
            dirInfo.D{IM}=inftMax;
        end
        
        IMold = IM;
    end    
    nTs(k)=IM;
    
end

impactTime = abs((caseInfo.dropCentre(2)-caseInfo.dropRadius)/caseInfo.dropVelocity(2));
[~,nT0]=min(abs(dirInfo.time-impactTime));
caseInfo.nT0=nT0;

dirInfo{:,varName}=cell(size(dirInfo,1),1);

t=caseInfo.nTi;
timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{t}];
[varC]=getVar(timePathName,dirInfo.geo{t}.centres,dirInfo.intf{t},varName,varDim);
dirInfo.(varName){t}=varC;
    
for k=1:nS
    t=nTs(k);   
    timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{t}];
    [varC]=getVar(timePathName,dirInfo.geo{t}.centres,dirInfo.intf{t},varName,varDim);
    dirInfo.(varName){t}=varC;
end

caseInfo.nTs=nTs;

end


function [caseInfo,dirInfo]=getDMax(caseInfo,dirInfo)

nT=nnz(~isnan([dirInfo.time]));

nTi=2;
nTmin=0;
nTmax = 0;
nTL = 0;
nTH = nT;
%nTM = floor((nTL+nTH)/2);

% Find Impact time
for t=2:nT %No data in t=1 or 0s
    timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{t}];
    geoInfo=getGeo(timePathName);
    dirInfo.geo{t}=geoInfo;
    [inftC,inftMax]=getIntfC(caseInfo,geoInfo,timePathName);
    dirInfo.intf{t}=inftC;
    dirInfo.D{t}=inftMax;
    if (inftMax>0)
        nTmin = t;
        nTL = nTmin;
        break
    end
end

if isempty(dirInfo.D{nTL})
    timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{nTL}];
    geoInfo=getGeo(timePathName);
    [inftC,inftMax]=getIntfC(caseInfo,geoInfo,timePathName);
    dirInfo.geo{nTL}=geoInfo;
    dirInfo.intf{nTL}=inftC;
    dirInfo.D{nTL}=inftMax;
end

if isempty(dirInfo.D{nTL+1})
    timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{nTL+1}];
    geoInfo=getGeo(timePathName);
    [inftC,inftMax]=getIntfC(caseInfo,geoInfo,timePathName);
    dirInfo.geo{nTL+1}=geoInfo;
    dirInfo.intf{nTL+1}=inftC;
    dirInfo.D{nTL+1}=inftMax;
end

if isempty(dirInfo.D{nTH})
    timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{nTH}];
    geoInfo=getGeo(timePathName);
    [inftC,inftMax]=getIntfC(caseInfo,geoInfo,timePathName);
    dirInfo.geo{nTH}=geoInfo;
    dirInfo.intf{nTH}=inftC;
    dirInfo.D{nTH}=inftMax;
end

if isempty(dirInfo.D{nTH-1})
    timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{nTH-1}];
    geoInfo=getGeo(timePathName);
    [inftC,inftMax]=getIntfC(caseInfo,geoInfo,timePathName);
    dirInfo.geo{nTH-1}=geoInfo;
    dirInfo.intf{nTH-1}=inftC;
    dirInfo.D{nTH-1}=inftMax;
end

if (dirInfo.D{nTL}>=dirInfo.D{nTL+1})
    nTmax = nTL;
elseif (dirInfo.D{nTH}>=dirInfo.D{nTH-1})
    nTmax = nTH;
else
    it = 0;
    while (nTL<=nTH)
        it =it+1;
      
        
        nTM = floor((nTL+nTH)/2);
        if isempty(dirInfo.D{nTM})
            timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{nTM}];
            geoInfo=getGeo(timePathName);
            [inftC,inftMax]=getIntfC(caseInfo,geoInfo,timePathName);
            dirInfo.geo{nTM}=geoInfo;
            dirInfo.intf{nTM}=inftC;
            dirInfo.D{nTM}=inftMax;
        end
        
        if isempty(dirInfo.D{nTM+1})
            timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{nTM+1}];
            geoInfo=getGeo(timePathName);
            [inftC,inftMax]=getIntfC(caseInfo,geoInfo,timePathName);
            dirInfo.geo{nTM+1}=geoInfo;
            dirInfo.intf{nTM+1}=inftC;
            dirInfo.D{nTM+1}=inftMax;
        end
        
        if isempty(dirInfo.D{nTM-1})
            timePathName = [caseInfo.fullPathName,filesep,dirInfo.name{nTM-1}];
            geoInfo=getGeo(timePathName);
            [inftC,inftMax]=getIntfC(caseInfo,geoInfo,timePathName);
            dirInfo.geo{nTM-1}=geoInfo;
            dirInfo.intf{nTM-1}=inftC;
            dirInfo.D{nTM-1}=inftMax;
        end
        
        if ((dirInfo.D{nTM}>=dirInfo.D{nTM+1}) &&  (dirInfo.D{nTM}>=dirInfo.D{nTM-1}))
            nTmax = nTM;
            break
        elseif (dirInfo.D{nTM}<dirInfo.D{nTM+1})
            nTL = nTM+1;
        else
            nTH = nTM-1;
        end
        
    end
end
caseInfo.nTi=nTi;
caseInfo.nTmin=nTmin;
caseInfo.nTmax=nTmax;
caseInfo.Dmax=dirInfo.D{nTmax};

end



%-------------------------

function [caseInfo]=getDropProperties(caseInfo)

filePath = [caseInfo.fullPathName,filesep,'constant',filesep,'transportProperties'];
fid=fopen(filePath,'rt');
tline=fgetl(fid);
tlines=cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline=fgetl(fid);
end
fclose(fid);

for n=1:size(tlines)
    if contains(tlines{n},{'dropRadius'})
        caseInfo.dropRadius = str2double(regexp(tlines{n}, '\d+\.\d+E[-+]?\d+', 'match'));
    end
    if contains(tlines{n},{'dropCentre'})
        newStr = erase(tlines{n},'dropCentre');
        newStr = erase(newStr,';');
        newStr = erase(newStr,'	');
        p=sort([strfind(newStr,'('),strfind(newStr,')'),strfind(newStr,' ')]);
        for m=1:length(p)-1
            v(m)=str2double(newStr(p(m)+1:p(m+1)-1));
        end
        caseInfo.dropCentre=v(~isnan(v));
    end
    if contains(tlines{n},{'dropVelocity'})
        newStr = erase(tlines{n},'dropVelocity');
        newStr = erase(newStr,';');
        newStr = erase(newStr,'	');
        p=sort([strfind(newStr,'('),strfind(newStr,')'),strfind(newStr,' ')]);
        for m=1:length(p)-1
            v(m)=str2double(newStr(p(m)+1:p(m+1)-1));
        end
        caseInfo.dropVelocity=v(~isnan(v));
    end
end

end

function [dirInfo]=findDirInfo(caseInfo)

dirInfo=dir(caseInfo.fullPathName);
if size(dirInfo,1)==0
    error('Empty folder')
end
dirInfo=dirInfo([dirInfo.isdir]==1);
dirInfo(ismember({dirInfo.name},{'.', '..'}))=[];

for k=1:size(dirInfo,1)
    dirInfo(k).time=str2num(dirInfo(k).name);
    if isempty(dirInfo(k).time)
        dirInfo(k).time=nan;
    end
end
[~,I]=sort([dirInfo.time]);

dirInfo=dirInfo(I,:);
dirInfo = struct2table(dirInfo);

end

function [geoInfo]=getGeo(timePathName)

geoPathName = [timePathName,filesep,'polyMesh'];

loopOut=cell(3,1);
loopVarName={'points','faces','owner'};
loopVarSize=[3,16,1];
parfor k=1:numel(loopOut)
    [tlines]=readOpenFoamFile([geoPathName,filesep,loopVarName{k}]);
    [val]=getValues(tlines,loopVarSize(k));
    loopOut{k}=val;
end
pts=loopOut{1};
faces=loopOut{2};
owner=loopOut{3};

C=unique(owner);
cells=cell(length(C),1);
for a = 1:length(C)
    idxF=find(a-1==owner);
    for b = 1:length(idxF)
        idxP=(faces(idxF(b),:)+1)';
        cells{a}=cat(1,cells{a},pts(idxP(~isnan(idxP)),:));
    end
end
geoInfo.cells=cells;

centres=nan(length(cells),2);
widths=nan(length(cells),2);
for a = 1:length(centres)
    centres(a,1)=(min(cells{a}(:,1))+max(cells{a}(:,1)))/2;
    centres(a,2)=(min(cells{a}(:,2))+max(cells{a}(:,2)))/2;
    widths(a,1)=max(cells{a}(:,1))-min(cells{a}(:,1));
    widths(a,2)=max(cells{a}(:,1))-min(cells{a}(:,1));
end

geoInfo.centres=centres;
geoInfo.widths=widths;

% colorP = hsv(length(C));
% for a = 1:length(C)
%     plot3(cells{a}(:,1),cells{a}(:,2),cells{a}(:,3),'.','color',colorP(a,:)), hold on
% end

end

function [varC]=getVar(timePathName,XY,inftC,varName,varDim)


filePath = [timePathName,filesep,varName];
[tlines]=readOpenFoamFile(filePath);
if isempty(tlines)
    
    valUni=zeros(1,varDim);
    
    fid=fopen(filePath,'rt');
    tline=fgetl(fid);
    tlines=cell(0,1);
    while ischar(tline)
        tlines{end+1,1} = tline;
        tline=fgetl(fid);
    end
    fclose(fid);
    
    for n=1:size(tlines)
        if contains(tlines{n},{'internalField'})
            if (varDim==1)
                valUni = str2double(regexp(tlines{n}, '\d+\.\d+', 'match'));
            else
                newStr = erase(tlines{n},'internalField');
                newStr = erase(newStr,'uniform');
                newStr = erase(newStr,';');
                newStr = erase(newStr,'	');
                p=sort([strfind(newStr,'('),strfind(newStr,')'),strfind(newStr,' ')]);
                for m=1:length(p)-1
                    v(m)=str2double(newStr(p(m)+1:p(m+1)-1));
                end
                valUni=v(~isnan(v));
            end
        end
    end
    
    varC=ones(size(inftC,1),1)*valUni;
    
else
    
    [varVal]=getValues(tlines,varDim);
    varC=nan(size(inftC,1),varDim);
    for k=1:size(inftC,1)
        dist2=sum((inftC(k,:) - XY).^2,2);
        [~,I]=sort(dist2,'ascend');
        I=I(1:2);
        coef1 = dist2(I(2))/(dist2(I(1))+dist2(I(2)));
        coef2 = dist2(I(1))/(dist2(I(1))+dist2(I(2)));
        
        varC(k,:)=coef1*varVal(I(1),:)+coef2*varVal(I(2),:);
        
    end
    
    varC(end,:)=varC(end-1,:);
end

end

function [inftC,inftMax]=getIntfC(caseInfo,geoInfo,timePathName)

[tlines]=readOpenFoamFile([timePathName,filesep,caseInfo.intfVar]);
[inftVal]=getValues(tlines,1);

tolX=max(unique(geoInfo.widths(:,1)))/10;
tolY=max(unique(geoInfo.widths(:,2)))/10;
dC = min(min(uniquetol(geoInfo.widths(:,1),tolX)),min(uniquetol(geoInfo.widths(:,1),tolY)));

%figure,plot3(geoInfo.centres(:,1),geoInfo.centres(:,2),inftVal,'*'), xticks(0:dC:4.5e-3), yticks(0:dC:4.5e-3), grid on, view(2), axis equal

%[X,Y] = meshgrid([dC/2:dC:max(geoInfo.centres(:,1))],[dC/2:dC:max(geoInfo.centres(:,2))]);
%dist2 = (X-reshape(geoInfo.centres(:,1),[1,1,length(geoInfo.centres(:,1))])).^2 + (Y-reshape(geoInfo.centres(:,2),[1,1,length(geoInfo.centres(:,2))])).^2;
%[~,idxMin] = min(dist2,[],3);
%Z = reshape(inftVal(idxMin),size(X));

% [X,Y] = meshgrid([dC/2:dC:max(geoInfo.centres(:,1))],[dC/2:dC:max(geoInfo.centres(:,2))]);
% distVal=pdist2([X(:) Y(:)], geoInfo.centres);
% [~, idxMin] = min(distVal, [], 2);
% Z = reshape(inftVal(idxMin),size(X));

[X,Y] = meshgrid([dC/2:dC:max(geoInfo.centres(:,1))],[dC/2:dC:max(geoInfo.centres(:,2))]);
Z = nan(size(X));

% subDomainNbr = caseInfo.subDomainNbr;
% num_elements = numel(Z);
% num_subgrids = round(num_elements / (num_elements/subDomainNbr)) - 1;

num_subgrids = caseInfo.subDomainNbr;

Z_sub = cell(num_subgrids, 1);
r_sub = cell(num_subgrids, 1);
c_sub = cell(num_subgrids, 1);

% num_rows = ceil(size(X, 1) / sqrt(num_subgrids));
% num_cols = ceil(size(X, 2) / sqrt(num_subgrids));
num_cols = floor(size(X, 2) / num_subgrids);

parfor i = 1:num_subgrids
    % Calculate the row and column indices for this sub-grid
%     row_start = floor((i-1) / sqrt(num_subgrids)) * num_rows + 1;
%     row_end = min(row_start + num_rows - 1, size(Z, 1));
%     col_start = floor(mod(i-1, sqrt(num_subgrids))) * num_cols + 1;
%     col_end = min(col_start + num_cols - 1, size(Z, 2));
    
%     if(row_end + num_rows/2 > size(Z, 1))
%         row_end = size(Z, 1);
%     end
%     if(col_end + num_cols/2 > size(X, 2))
%         col_end = size(Z, 2);
%     end
    
%     adjusted_rows = row_start:row_end;
%     adjusted_cols = col_start:col_end;
%     r_sub{i} = adjusted_rows;
%     c_sub{i} = adjusted_cols;
    

    adjusted_rows = 1:size(Z, 1);
    col_start = (i-1) * num_cols + 1;
    col_end = min(i * num_cols, size(X, 2));
    if(col_end + num_cols/2 > size(X, 2))
        col_end = size(Z, 2);
    end
    adjusted_cols = col_start:col_end;
    r_sub{i} = adjusted_rows;
    c_sub{i} = adjusted_cols;
    
    % Extract the sub-grid
    subgrid_X = X(adjusted_rows, adjusted_cols);
    subgrid_Y = Y(adjusted_rows, adjusted_cols);
    
    gridMargin = dC*10;
    idxSubX = ((geoInfo.centres(:,1)>min(subgrid_X(:))-gridMargin) .* (geoInfo.centres(:,1)<max(subgrid_X(:))+gridMargin) >0.5);
    idxSubY = ((geoInfo.centres(:,2)>min(subgrid_Y(:))-gridMargin) .* (geoInfo.centres(:,2)<max(subgrid_Y(:))+gridMargin) >0.5);
    idxSubXY = (idxSubX .* idxSubY >0.5);
    
    subCentres = geoInfo.centres(idxSubXY,:);
    subInftVal = inftVal(idxSubXY);
    
    %Associate values
    distVal = pdist2([subgrid_X(:) subgrid_Y(:)], subCentres); %distVal = pdist2([subgrid_X(:) subgrid_Y(:)], geoInfo.centres);
    
    [~, idxMin] = min(distVal, [], 2);
    
    sub_Z = nan(length(adjusted_rows), length(adjusted_cols));
    sub_Z(:) = subInftVal(idxMin); %sub_Z(:) = inftVal(idxMin);
    Z_sub{i} = sub_Z;
end

for i = 1:num_subgrids
    Z(r_sub{i},c_sub{i}) = Z_sub{i};
end

XMod = cat(1,X(1,:),X);
YMod = cat(1,-Y(1,:),Y);
ZMod = cat(1,Z(1,:),Z);
XMod = cat(2,-XMod(:,1),XMod);
YMod = cat(2,YMod(:,1),YMod);
ZMod = cat(2,ZMod(:,1),ZMod);
%figure,plot3(geoInfo.centres(:,1),geoInfo.centres(:,2),inftVal,'*'), hold on, scatter3(XMod(:),YMod(:),ZMod(:),5,ZMod(:)),view(2), axis equal
%figure,plot3(geoInfo.centres(:,1),geoInfo.centres(:,2),inftVal,'*'), hold on, surf(XMod,YMod,ZMod),view(2), axis equal

C=contourc(XMod(1,:),YMod(:,1),ZMod,[caseInfo.intfLim caseInfo.intfLim])'; %C=contour(XMod,YMod,ZMod,[0.5 0.5])';

%C=contourc(X(1,:),Y(:,1),Z,[0.5 0.5])'; %C=contour(X,Y,Z,[0.5 0.5])';

%Select longest section
i=1;
while i<size(C,1)
    iNew=i+C(i,2)+1;
    C(i,:)=nan(1,2);
    i=iNew;
end
C=cat(1,C,nan(1,2));
nanPos=find(isnan(C(:,1)));
[~,idx]=max(diff(nanPos));
C=C(nanPos(idx)+1:nanPos(idx+1)-1,:);

%Contour orientation 
flipOn=false;
if (C(1,1)*C(end,1)>0)
    if (C(1,2)<C(end,2))
        flipOn = true;
    end
else
    if (C(1,1)>C(end,1))
        flipOn = true;
    end
end
if (flipOn==true)
	C=flipud(C);
end

% Extrapolation

if (C(1,1)<0)
    C(1,1)=0;
    C(1,2)=C(2,2);
end
if (C(end,1)<0)
    C(end,1)=0;
    C(end,2)=C(end-1,2);
end

if (C(1,2)<0)
    C(1,1)=C(2,2);
    C(1,2)=0;
end
if (C(end,2)<0)
    C(end,1)=C(end-1,1);
    C(end,2)=0;
    
    try
        [varC]=getVar(timePathName,geoInfo.centres,C,caseInfo.CAVar,1);
        m = tand(180-varC(end-1));
        p = C(end-1,2)-m*C(end-1,1);
        x0 = -p/m;
        if x0<0
            x0=0;
        end
        C(end,1)=x0;
    catch
        %disp('Contact angle ignored')
    end
    
end

inftC = C;
inftMax = C(end,1);

end



function [tlines]=readOpenFoamFile(filePath)
fid=fopen(filePath,'rt');
tline=fgetl(fid);
tlines=cell(0,1);
while ischar(tline)
    tlines{end+1,1} = tline;
    tline=fgetl(fid);
end
fclose(fid);
tlines=tlines(find(ismember(tlines,{'('})==1,1,'first')+1:1:find(ismember(tlines,{')'})==1,1,'first')-1);
end

function [vals]=getValues(tlines,dimVal)
vals=nan(length(tlines),dimVal);
if dimVal==1
    for n=1:length(tlines)
        vals(n)=str2double(tlines{n});
    end
else
    for n=1:length(tlines)
        p=sort([strfind(tlines{n},'('),strfind(tlines{n},')'),strfind(tlines{n},' ')]);
        for m=1:length(p)-1
            vals(n,m)=str2double(tlines{n}(p(m)+1:p(m+1)-1));
        end
    end
end
end
