% ring connection statistics
%% Distribution of Shared Ring Number for Each Neuron
% if there should be some limitation condition concerning the ring showup
% times, or there is no exhaustion method (except a exhaustion set of
% neurons)
% standard : at least 1 showup in one trial in average
% temp standard: >10 showup in the total recording session
load rings.mat
load ringSampleTrials.mat

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

sessList = [metadata{1}(:,1);metadata{2}(:,1);metadata{3}(:,1)];
sessList = unique(sessList);

%% save each session and sample trials
sampleTrials = cell(1,3);
sampleTrials{3} = s1s2Trials;
save('ringSampleTrials.mat','sampleTrials');

%%
% su1 ... su5(if possible)
ringConnListFull = cell(114,1);
 
tmp = [-1,-1,-1,-1,-1];
for midx = 1:3
    load(sprintf('turnStatOldLen%dFullSample.mat',midx+2));
    turnStat = turnStatOld;
    %load(sprintf('turnStatLen%dFullSample.mat',midx+2));
    %load(sprintf('turnStatLen%dFull.mat',midx+2));
for i = sessList(:)'
    trials = [sampleTrials{1}(metadata{1}(:,1) == i,:);sampleTrials{2}(metadata{2}(:,1) == i,:);sampleTrials{3}(metadata{3}(:,1) == i,:)];
    if isempty(samples{midx}(metadata{midx}(:,1)==i))
        continue;
    end
    currSamples = mode(samples{midx}(metadata{midx}(:,1)==i));
    trials1 = trials{1,1};
    trials2 = trials{1,2};
    tnum = length(trials1)+length(trials2);
    %if currSamples == 1
    %    tnum = length(trials{1,1});
    %elseif currSamples == 2
    %    tnum = length(trials{1,2});
    %else
    %    disp('check!')
    %    pause;
    %end
    currRingsStat = turnStat(metadata{midx}(:,1)==i);
    currRings = metadata{midx}(metadata{midx}(:,1)==i,2:end);
    for j = 1:size(currRings,1)
        counts = sum(cell2mat(currRingsStat{j}'));
        if counts >= tnum
            currRingTemp = tmp;
            currRingTemp(1:midx+2) = currRings(j,:);
            ringConnListFull{i} = [ringConnListFull{i};currRingTemp];
        end
    end
end
end

% shuffle data



sharedRingRate = [];
for ssidx = 1:length(ringConnListFull)
    tempRings = ringConnListFull{ssidx};
    if size(tempRings,1) < 10
        continue;
    end
    suList = unique(tempRings);
    suList(suList==-1) = [];
    for suidx = 1:length(suList)
        shareRingCounts = length(find(tempRings==suList(suidx)));
        sharedRingRate = [sharedRingRate;shareRingCounts/size(tempRings,1)];
    end
end
hist(sharedRingRate,200)
ylabel('Probability');
xlabel('Fraction of rings the neuron take part in');

    
    

%% Distribution of Ring Number Sharing Neuron for Each Ring
shareSuRingRate = [];
for ssidx = 1:length(ringConnListFull)
    tempRings = ringConnListFull{ssidx};
    if size(tempRings,1)< 10
        continue;
    end
    for ridx = 1:size(tempRings,1)
        shareSuRingCount = sum(any(tempRings == tempRings(ridx,:) & tempRings ~= -1,2))-1;
        shareSuRingRate = [shareSuRingRate;shareSuRingCount/size(tempRings,1)];
    end
end
hist(shareSuRingRate,200);
ylabel('Probability');
xlabel('Fraction of rings sharing at least one same neuron');


%% Fraction of Rings have shared neuron in 1s bin

    