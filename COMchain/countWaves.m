% function countWaves(opt)
% % prepare functional coupling linked neuronal-wave data
% arguments
%     opt.peak (1,1) logical = false
%     opt.strict (1,1) logical = true %strict FC criteria
% %     opt.min_diff (1,1)
%     opt.screen (1,1) logical = false
% end
opt.peak = false;
opt.strict =true;
chainsumCE = zeros(116,4,4); % S1 C S2 C S1 E S2E (2,3,4,5)
realChainsCE = zeros(116,4,4); % S1 C S2 C S1 E S2E (2,3,4,5)
realChainsCELabel = cell(116,4,4); % S1 C S2 C S1 E S2E (2,3,4,5)
% rate_congru = 0.018;
% rate_incongru = 0.015;
% rate_nonmem = 0.012;
load('../sums_conn.mat','sums_conn_str');
chainsum = cell(length(sums_conn_str),2); % correct error
% chainsum = cell(length(sums_conn_str),6); % S1 congru, S2 congru, incongru, non-mem, S1 congru error, S2 congru error
% chainsumPred = cell(length(sums_conn_str),6); % S1 congru, S2 congru, incongru, non-mem, S1 congru error, S2 congru error
for chainlen = 2:5
for fidx=1:116
    disp(fidx)
    if opt.strict
        ccgqc=sums_conn_str(fidx).qc; %reference quality control parameter
        strict_sel=ccgqc(:,1)>0; %& ccg(:,4)<350 & ccg(:,4)>252 & ccg(:,6)<20 & ccg(:,5)<20 & ccg(:,6)<20 ;
        %1:Polar 2:Time of peak 3:Noise peaks 4:FWHM
        oneccg=sums_conn_str(fidx).ccg_sc(strict_sel,:); %ccg
        onecon=sums_conn_str(fidx).sig_con(strict_sel,:); %jitter controlled significant functional coupling    
    else % full input from English, Buzsaki code
        oneccg=sums_conn_str(fidx).ccg_sc;
        onecon=sums_conn_str(fidx).sig_con;
    end
    %onecom=wave.get_com_map('onepath',strrep(sums_conn_str(fidx).folder,'zx/neupix','xd/data'),'peak',opt.peak); %per su center of mass, normalized FR
    load(sprintf('com_str_%d.mat',fidx));
    % construct data for sorting by com
    % compare correct error trials
    skey=fieldnames(com_str);
%     nonsel=size(onecom.(skey{1}).s0,1);
%     sel=length(onecom.(skey{1}).s1) + length(onecom.(skey{1}).s2);
%     selmax = [max([nonsel,selmax(1)]),max([sel,selmax(2)])];
%     disp([fidx,nonsel,sel])
    % exist pair        
    if ~isempty(skey) && ~isempty(com_str.(skey{1}).s1e) && ~isempty(com_str.(skey{1}).s2e)
        % correct 
        s1sel=cell2mat(com_str.(skey{1}).s1.keys); %pre-selected transient selective su
        s2sel=cell2mat(com_str.(skey{1}).s2.keys);
        s1com=cell2mat(com_str.(skey{1}).s1.values);
        s2com=cell2mat(com_str.(skey{1}).s2.values);
        preSortMat=[s1com,s2com;ones(1,length(s1com)),2*ones(1,length(s2com));double(s1sel),double(s2sel)]';
        sortedMat=sortrows(preSortMat);
        s1cdata=sortedMat(sortedMat(:,2)==1,:);
        s1clabel=uint16(s1cdata(:,3));
        s2cdata=sortedMat(sortedMat(:,2)==2,:);
        s2clabel=uint16(s2cdata(:,3));
        if length(s1clabel)<=intmax('uint8')
            s1Idx = uint8(1:length(s1clabel));
        else
            s1Idx = uint16(1:length(s1clabel));
        end
        if length(s1clabel)<=intmax('uint8')
            s2Idx = uint8(1:length(s2clabel));
        else
            s2Idx = uint16(1:length(s2clabel));
        end
        
        idxperms = nchoosek(s1Idx,chainlen);
        timeperms = arrayfun(@(x) sortedMat(x,1),idxperms);
        labelperms = arrayfun(@(x) s1clabel(x),idxperms);
        clear idxperms
        timediff = diff(timeperms,1,2);
        clear timeperms
        timesel = all(timediff >50/250 & timediff<=8,2);
        clear timediff
        allperms = labelperms(timesel,:);
        clear timesel labelperms
        chainsumCE(fidx,1,chainlen-1) = size(allperms,1);
        selector = ones(size(allperms,1),1);
        for level = 1:chainlen-1
            selector = selector & ismember(allperms(:,level:level+1),onecon,'rows');
        end
        realChainsCELabel{fidx,1,chainlen-1} = allperms(selector,:);
        clear allperms
        realChainsCE(fidx,1,chainlen-1) = size(realChainsCELabel{fidx,1,chainlen-1},1);
        
        idxperms = nchoosek(s2Idx,chainlen);
        timeperms = arrayfun(@(x) sortedMat(x,1),idxperms);
        labelperms = arrayfun(@(x) s2clabel(x),idxperms);
        clear idxperms
        timediff = diff(timeperms,1,2);
        clear timeperms
        timesel = all(timediff >50/250 & timediff<=8,2);
        clear timediff
        allperms = labelperms(timesel,:);
        clear timesel labelperms
        chainsumCE(fidx,2,chainlen-1) = size(allperms,1);
        selector = ones(size(allperms,1),1);
        for level = 1:chainlen-1
            selector = selector & ismember(allperms(:,level:level+1),onecon,'rows');
        end
        realChainsCELabel{fidx,2,chainlen-1} = allperms(selector,:);
        clear allperms
        realChainsCE(fidx,2,chainlen-1) = size(realChainsCELabel{fidx,2,chainlen-1},1);
        % error
        s1sel=cell2mat(com_str.(skey{1}).s1e.keys); %pre-selected transient selective su
        s2sel=cell2mat(com_str.(skey{1}).s2e.keys);
        s1com=cell2mat(com_str.(skey{1}).s1e.values);
        s2com=cell2mat(com_str.(skey{1}).s2e.values);
        preSortMat=[s1com,s2com;ones(1,length(s1com)),2*ones(1,length(s2com));double(s1sel),double(s2sel)]';
        sortedMat=sortrows(preSortMat);
        s1cdata=sortedMat(sortedMat(:,2)==1,:);
        s1clabel=uint16(s1cdata(:,3));
        s2cdata=sortedMat(sortedMat(:,2)==2,:);
        s2clabel=uint16(s2cdata(:,3));
        if length(s1clabel)<=intmax('uint8')
            s1Idx = uint8(1:length(s1clabel));
        else
            s1Idx = uint16(1:length(s1clabel));
        end
        if length(s1clabel)<=intmax('uint8')
            s2Idx = uint8(1:length(s2clabel));
        else
            s2Idx = uint16(1:length(s2clabel));
        end
        
        idxperms = nchoosek(s1Idx,chainlen);
        timeperms = arrayfun(@(x) sortedMat(x,1),idxperms);
        labelperms = arrayfun(@(x) s1clabel(x),idxperms);
        clear idxperms
        timediff = diff(timeperms,1,2);
        clear timeperms
        timesel = all(timediff >50/250 & timediff<=8,2);
        clear timediff
        allperms = labelperms(timesel,:);
        clear timesel labelperms
        chainsumCE(fidx,3,chainlen-1) = size(allperms,1);
        selector = ones(size(allperms,1),1);
        for level = 1:chainlen-1
            selector = selector & ismember(allperms(:,level:level+1),onecon,'rows');
        end
        realChainsCELabel{fidx,3,chainlen-1} = allperms(selector,:);
        clear allperms
        realChainsCE(fidx,3,chainlen-1) = size(realChainsCELabel{fidx,3,chainlen-1},1);
        
        idxperms = nchoosek(s2Idx,chainlen);
        timeperms = arrayfun(@(x) sortedMat(x,1),idxperms);
        labelperms = arrayfun(@(x) s2clabel(x),idxperms);
        clear idxperms
        timediff = diff(timeperms,1,2);
        clear timeperms
        timesel = all(timediff >50/250 & timediff<=8,2);
        clear timediff
        allperms = labelperms(timesel,:);
        clear timesel labelperms
        chainsumCE(fidx,4,chainlen-1) = size(allperms,1);
        selector = ones(size(allperms,1),1);
        for level = 1:chainlen-1
            selector = selector & ismember(allperms(:,level:level+1),onecon,'rows');
        end
        realChainsCELabel{fidx,4,chainlen-1} = allperms(selector,:);
        clear allperms
        realChainsCE(fidx,4,chainlen-1) = size(realChainsCELabel{fidx,4,chainlen-1},1);
 
    end
    save(sprintf('chainSumAll_%d_%d.mat',fidx,chainlen),'chainsumCE','realChainsCE','realChainsCELabel','-v7.3');
end
end

%     for i = 2:length(realChains)
%         currChain = ID{i};
%         selector = ones(size(currChain,1),1);
%         for level = 1:size(currChain,2)-1
%             selector = selector & ismember(currChain(:,level:level+1),onecon,'rows');
%         end
%         realChain = currChain(selector,:);
        
 
    

        
        
% function [repertoire,sel,ID] = permschain(sortedMat,chainlen)
% tic
% p = 0.02; % Buzaki maximum ratio
% repertoire = cell(1,chainlen);
% sel = cell(1,chainlen);
% tmax = sortedMat(end,1);
% depth = 1;
% repertoire{1} = sortedMat(1:size(sortedMat,1)-chainlen,1);
% sel{1} = sortedMat(1:size(sortedMat,1)-chainlen,2);
% ID{1} = sortedMat(1:size(sortedMat,1)-chainlen,3);
% loopIdx = 1;
% while depth < chainlen
%     while loopIdx <= size(repertoire{depth},1)
%         tpre = repertoire{depth}(loopIdx,:);
%         selpre = sel{depth}(loopIdx,:);
%         IDpre = ID{depth}(loopIdx,:);
%         tail = sortedMat(sortedMat(:,1)>tpre(end)+50/250 & sortedMat(:,1)<=tpre(end)+8 & sortedMat(:,1)<=tmax,:);
%         repertoire{depth+1} = [repertoire{depth+1};ones(length(tail(:,1)),1)*tpre,tail(:,1)];
%         sel{depth+1} = [sel{depth+1};ones(length(tail(:,2)),1)*selpre,tail(:,2)];
%         ID{depth+1} = [ID{depth+1};ones(length(tail(:,3)),1)*IDpre,tail(:,3)];
%         loopIdx = loopIdx + 1;
%     end
%     if isempty(repertoire{depth+1})
%         break;
%     end
%     if size(repertoire{depth+1},1)*(p^depth)<1
%         break;
%     end
%     depth = depth + 1;
%     loopIdx = 1;
% end
% toc
% end


% function [] = identicalSelChains()
% selorder = [1,1,2,1,1];
% p = 0.02;
% repertoire = cell(1,length(selorder));
% loopIdx = 1;
% tmax = sortedMat(end,1);
% depth = 1;
% repertoire{1} = sortedMat(sortedMat(:,2)==selorder(1),1);
% while depth < length(selorder)
%     while loopIdx <= size(repertoire{depth},1)
%         tpre = repertoire{depth}(loopIdx,:);
%         tail = sortedMat(sortedMat(:,1)>tpre(end)+50/250 & sortedMat(:,1)<=tpre(end)+8 & sortedMat(:,1)<=tmax & sortedMat(:,2)==selorder(depth+1),:);
%         repertoire{depth+1} = [repertoire{depth+1};ones(length(tail(:,1)),1)*tpre,tail(:,1)];
%         loopIdx = loopIdx + 1;
%     end
%     if isempty(repertoire{depth+1})
%         break;
%     end
%     if size(repertoire{depth+1},1)*(p^depth)<1
%         break;
%     end
%     depth = depth + 1;
%     loopIdx = 1;
% end
        



