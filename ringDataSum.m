% rings data sum for loose criteria
%idces from hebb_pattern_showcase.m
if ~exist('rings','var')
    load rings.mat
    load 114_sorted_file_path.mat
    %freg=load('sel_conn_chain_duo_6s_1_2.mat','pair_reg','pair_chain');
    %load reg_keep.mat
    cd('/home/xd/pkg')
    addpath(fullfile('npy-matlab-master','npy-matlab'));
    addpath('fieldtrip-20200521');
    ft_defaults
    delay=6;
    if isunix
        homedir='/home/zx/neupix/wyt';
    elseif ispc
        homedir='k:\neupix\wyt';
    end
end

ringsm = cell(1,3);
metadata = cell(1,3);
samples = cell(1,3);
for midx=1:3
r1=[];
for bb=1:6
    r1=[r1;cell2mat(rings(midx,:,bb,1)')];
end
r2=[];
for bb=1:6
    r2=[r2;cell2mat(rings(midx,:,bb,2)')];
end
r1(:,midx+3)=1;r2(:,midx+3)=2;ringsm{midx}=[r1(:,1:(midx+2));r2(:,1:(midx+2))];samps=[r1(:,midx+3);r2(:,midx+3)];
samples{midx} = samps;
metadata{midx} = zeros(length(ringsm{midx}),midx+3);
for ssidx=1:length(ringsm{midx})
    sessIdx=idivide(int32(ringsm{midx}(ssidx,1)),int32(100000));
    transIds=int32(rem(ringsm{midx}(ssidx,:),100000))';
    metadata{midx}(ssidx,1) = sessIdx;
    metadata{midx}(ssidx,2:end) = transIds';
end
end
ringList = -1*ones(length(metadata{1})+length(metadata{2})+length(metadata{3}),6);
ringList(1:length(metadata{1}),1:4) = metadata{1};
ringList(length(metadata{1})+1:length(metadata{1})+length(metadata{2}),1:5) = metadata{2};
ringList(length(metadata{1})+length(metadata{2})+1:length(metadata{1})+length(metadata{2})+length(metadata{3}),1:6) = metadata{3};
sessList = [metadata{1}(:,1);metadata{2}(:,1);metadata{3}(:,1)];
sessList = unique(sessList);


ringsDset = cell(114,1);
for i = sessList(:)'
    folder=regexp(sorted_fpath{i},'(\w|\\|-)*','match','once');
    [folderType,file,spkFolder,metaFolder,~]=jointFolder(folder,cell(0),homedir);
    currmodel='selec';
    sustIds=[];nonselIds=[];
    currRings = ringList(ringList(:,1)==i,2:end);
    ringsDset{i} = [];
    for ridx = 1:size(currRings,1)
        currRing = currRings(ridx,:);
        msize = length(find(currRing ~= -1));
        transIds = currRing(1:msize)';
        [avail,spktrial]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,currmodel);
        cfg.trials = find(spktrial.trialinfo(:,8)==delay);
        for tidx=cfg.trials(:)'
            ts_id=[];
            for seqid=1:msize
                tsel=spktrial.trial{seqid}==tidx;
                ts=spktrial.time{seqid}(tsel);
                ts=ts(ts>1 & ts<7);
                ts(2,:)=seqid;
                ts_id=[ts_id;ts'];
            end
            if length(ts_id) < msize + 1
                continue;
            end
            [~,s]=sort(ts_id(:,1));
            ts_id=ts_id(s,:);
            tempdatasum=relax_tag(ts_id,msize);
            if isempty(tempdatasum)
                continue;
            end
            ringNum = size(tempdatasum,1);
            tempInfoMat = -1*ones(ringNum,11);
            tempInfoMat(:,1:4) = tempdatasum;
            tempInfoMat(:,5) = tidx;
            tempInfoMat(:,6) = msize;
            tempInfoMat(:,7:6+msize) = ones(ringNum,1)*transIds(:)';
            ringsDset{i} = [ringsDset{i};tempInfoMat];
        end
    end
    if mod(i,5) == 0
        save('ringsDataSum.mat','ringsDset','i','-v7.3');
    end
end
save('ringsDataSum.mat','ringsDset','i','-v7.3');

function [folderType,file,spkFolder,metaFolder,error_list]=jointFolder(folder,error_list,homedir)
metaFolder=replace(folder,'\','/');
metaFolder=fullfile(fullfile(homedir,'DataSum'),metaFolder);
if isfolder(metaFolder)
    spkFolder=replace(metaFolder,'imec1','imec0');
    file=dir(fullfile(spkFolder,'spike_info.mat'));
    if isempty(file)
        folderType=-1;
        file=[];
        spkFolder=[];
        disp('Error processing file 2-tracks');
        disp(metaFolder);
        error_list(end+1,:)={folderType,metaFolder};
        %             pause;
        return
    end
    folderType=2;
else
    metaFolder=replace(metaFolder,'DataSum','DataSum/singleProbe');
    spkFolder=metaFolder;
    file=dir(fullfile(spkFolder,'spike_times.npy'));
    if isempty(file)
        folderType=-1;
        file=[];
        spkFolder=[];
        disp('Error processing file 1-track');
        disp(metaFolder);
        error_list(end+1,:)={folderType,metaFolder};
        return
    end
    folderType=1;
end
end
    
    
function [avail,out]=pre_process(folderType,spkFolder,metaFolder,sustIds,transIds,nonselIds,model)
sps=30000;
trials=clearBadPerf(h5read(fullfile(metaFolder,'events.hdf5'),'/trials')',model);
if isempty(trials)
    avail=false;
    out=[];
    return
end

%     trials=double(trials);
%     info=[trials(:,1)/s1s,trials(:,2)/s1s,trials(:,5),trials(:,6),trials(:,7),trials(:,8)];
if strcmp(model, 'full')
    cluster_ids=[sustIds;transIds;nonselIds];
elseif startsWith(model,'selec')
    cluster_ids=[sustIds;transIds];
elseif startsWith(model,'nonsel')
    cluster_ids=nonselIds;
else
    keyboard
end

%  single-unit candidate

if folderType==1
    spkTS=readNPY(fullfile(spkFolder,'spike_times.npy'));
    spkId=readNPY(fullfile(spkFolder,'spike_clusters.npy'));
elseif folderType==2
    fstr=load(fullfile(spkFolder,'spike_info.mat'));
    spkId=double([fstr.spike_info{1}{1};fstr.spike_info{1}{2}]);
    spkTS=double([fstr.spike_info{2}{1};fstr.spike_info{2}{2}]);
end
FT_SPIKE=struct();

FT_SPIKE.label=strtrim(cellstr(num2str(cluster_ids)));
FT_SPIKE.timestamp=cell(1,numel(cluster_ids));
for i=1:numel(cluster_ids)
    FT_SPIKE.timestamp{i}=spkTS(spkId==cluster_ids(i))';
end
%  continuous format F T struct file
cfg=struct();
cfg.trl=[trials(:,1)-3*sps,trials(:,1)+11*sps,zeros(size(trials,1),1)-3*sps,trials];
cfg.trlunit='timestamps';
cfg.timestampspersecond=sps;

FT_SPIKE=ft_spike_maketrials(cfg,FT_SPIKE);

out=FT_SPIKE;
avail=true;
end



function out=clearBadPerf(facSeq, mode)
if strcmp(mode, 'error')
    if length(facSeq)>=40
        errorsel=~xor(facSeq(:,5)==facSeq(:,6) , facSeq(:,7)>0);
        out=facSeq(errorsel,:);
    else
        out=[];
    end
else
    if length(facSeq)>=40
        facSeq(:,9)=0;
        i=40;
        while i<=length(facSeq)
            good=xor(facSeq(i-39:i,5)==facSeq(i-39:i,6) , facSeq(i-39:i,7)>0);
            facSeq(i-39:i,10)=good;
            if nnz(good)>=30 %.75 correct rate
                facSeq(i-39:i,9)=1;
            end
            i=i+1;
        end
        out=facSeq(all(facSeq(:,9:10),2),:);
    else
        out=[];
    end
end
end


function out=relax_tag(in,msize)
out = [];
skiptag=0;
for i=1:length(in)
    if i<skiptag
        continue
    end
    curr_su=in(i,2);
    targets=[(curr_su+1):(curr_su+msize-1),curr_su];
    targets(targets>msize)=targets(targets>msize)-msize;
    tsseq=[i,in(i,1:2)];
    for t=targets
        rows=tsseq(end,1)+(1:msize*10);
        rows(rows>length(in))=[];
        if isempty(rows)
            break
        end
        didx=find( ...
            in(rows,2)==t ... %post unit
            & in(rows,1)<tsseq(end,2)+0.01 ...
            & in(rows,1)>tsseq(end,2)+0.0005 ... %matching time window, assuming 1kHz
            ,1);
        if isempty(didx)
            break
        else
            tsseq=[tsseq;tsseq(end,1)+didx,in(tsseq(end,1)+didx,1:2)];
        end
    end
    if length(tsseq)<msize+1
        continue
    else
        out = [out;in(i,2),in(i,1),tsseq(end,3),tsseq(end,2)];
        skiptag=tsseq(2,1);
    end
end
end