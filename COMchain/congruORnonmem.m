% % function countWaves(opt)
% % prepare functional coupling linked neuronal-wave data
% arguments
%     opt.peak (1,1) logical = false
%     opt.strict (1,1) logical = true %strict FC criteria
% %     opt.min_diff (1,1)
%     opt.screen (1,1) logical = false
% end
opt.peak = false;
opt.strict =true;
chainsumCE = zeros(116,2,4); % nonmem incongru (2,3,4,5)
realChainsCE = zeros(116,2,4); % nonmem incongru (2,3,4,5)
realChainsCELabel = cell(116,2,4); % nonmem incongru (2,3,4,5)
% rate_congru = 0.018;
% rate_incongru = 0.015;
% rate_nonmem = 0.012;
load('../sums_conn.mat','sums_conn_str');
chainsum = cell(length(sums_conn_str),2); % correct error
% chainsum = cell(length(sums_conn_str),6); % S1 congru, S2 congru, incongru, non-mem, S1 congru error, S2 congru error
% chainsumPred = cell(length(sums_conn_str),6); % S1 congru, S2 congru, incongru, non-mem, S1 congru error, S2 congru error
for chainlen = 2:3
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
    load(sprintf('com_str_Newcrit_%d.mat',fidx));
    onecom = com_str;
    % construct data for sorting by com
    % compare correct error trials
    skey=fieldnames(onecom);
%     nonsel=size(onecom.(skey{1}).s0,1);
%     sel=length(onecom.(skey{1}).s1) + length(onecom.(skey{1}).s2);
%     selmax = [max([nonsel,selmax(1)]),max([sel,selmax(2)])];
%     disp([fidx,nonsel,sel])
    % exist pair 
    if ~isempty(skey)
        nonsel=cell2mat(com_str.(skey{1}).s0.keys);
        s1sel=cell2mat(com_str.(skey{1}).s1.keys); %pre-selected transient selective su
        s2sel=cell2mat(com_str.(skey{1}).s2.keys);
        s0com=cell2mat(com_str.(skey{1}).s0.values);
        s1com=cell2mat(com_str.(skey{1}).s1.values);
        s2com=cell2mat(com_str.(skey{1}).s2.values);
        nonmemPreSortMat = [s0com;double(nonsel)]';
        sortedNonmem = sortrows(nonmemPreSortMat);
        nonmemLabel=uint16(nonsel);
        preSortMat=[s1com,s2com;ones(1,length(s1com)),2*ones(1,length(s2com));double(s1sel),double(s2sel)]';
        sortedMat=sortrows(preSortMat);
        incongruLabel=uint16(sortedMat(:,3));
        if length(nonmemLabel)<=intmax('uint8')
            nonselIdx = uint8(1:length(nonmemLabel));
        else
            nonselIdx = uint16(1:length(nonmemLabel));
        end
        if length(incongruLabel)<=intmax('uint8')
            incongruIdx = uint8(1:length(incongruLabel));
        else
            incongruIdx = uint16(1:length(incongruLabel));
        end
        % non-mem
        idxperms = nchoosek(nonselIdx,chainlen);
        timeperms = arrayfun(@(x) sortedNonmem(x,1),idxperms);
        labelperms = arrayfun(@(x) nonmemLabel(x),idxperms);
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
        % incongru
        idxperms = nchoosek(incongruIdx,chainlen);
        typeperms = arrayfun(@(x) sortedMat(x,2),idxperms);
        incongrusel = any(diff(typeperms,1,2),2);
        clear typeperms incongrusel
        idxperms = idxperms(incongrusel);
        timeperms = arrayfun(@(x) sortedMat(x,1),idxperms);
        labelperms = arrayfun(@(x) incongruLabel(x),idxperms);
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
    end

save(sprintf('chainSumAllNcrit_%d_%d.mat',fidx,chainlen),'chainsumCE','realChainsCE','realChainsCELabel','-v7.3');

    
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
        



