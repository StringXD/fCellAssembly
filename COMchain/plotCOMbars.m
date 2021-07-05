
% compare chains between correct error
ratio = realChainsCE./chainsumCE;
conn = realChainsCE(:,:,2);
pairs = chainsumCE(:,:,2);
% correct error
a = (conn(:,1)+conn(:,2))./(pairs(:,1)+pairs(:,2)) ;
b = (conn(:,3)+conn(:,4))./(pairs(:,3)+pairs(:,4));
% nonmem incongru
%a = conn(:,1)./pairs(:,1) ;
%b = conn(:,2)./pairs(:,2);

sel = all(pairs>100,2);

meanCOMconnC = mean(a(sel));
meanCOMconnE = mean(b(sel));

semCOMconnC = std(a(sel))/sqrt(sum(sel));
semCOMconnE = std(b(sel))/sqrt(sum(sel));

fh=figure('Color','w','Position',[100,100,235,270]);
hold on
bar([meanCOMconnC,meanCOMconnE,meanCOMNM,meanCOMIncongru]*100,'FaceColor','w')
errorbar(1:4,[meanCOMconnC,meanCOMconnE,meanCOMNM,meanCOMIncongru]*100,[semCOMconnC,semCOMconnE,semCOMNM,semCOMIncongru]*100,'k.','CapSize',12)
ax = gca;
%ax.YLim = [0,1.5];
ax.YLim = [0, 0.05];
set(gca,'XTick',1:4,'XTickLabel',{'Correct','Error','Nonmem','Incongru'},'XTickLabelRotation',30)
ylabel('% 3 Neuron Wave')

%% diveded by FC pairs
% 2 neuron wave
% typenum = zeros(116,4);
% sumcom = zeros(116,4);
% for fidx = 1:116
%     ccgqc=sums_conn_str(fidx).qc; %reference quality control parameter
%     strict_sel=ccgqc(:,1)>0; %& ccg(:,4)<350 & ccg(:,4)>252 & ccg(:,6)<20 & ccg(:,5)<20 & ccg(:,6)<20 ;
%         %1:Polar 2:Time of peak 3:Noise peaks 4:FWHM
%     oneccg=sums_conn_str(fidx).ccg_sc(strict_sel,:); %ccg
%     onecon=sums_conn_str(fidx).sig_con(strict_sel,:); 
%     load(sprintf('com_str_%d.mat',fidx));
%     skey=fieldnames(com_str);
%     if ~isempty(skey) && ~isempty(com_str.(skey{1}).s1e) && ~isempty(com_str.(skey{1}).s2e)
%         for type = 1:4
%             switch(type)
%                 case 1
%                     samp = "s1";
%                 case 2
%                     samp = "s2";
%                 case 3
%                     samp = "s1e";
%                 case 4
%                     samp = "s2e";
%             end
%             mapkeys=cell2mat(com_str.(skey{1}).(samp).keys);
%             typesel=all(ismember(int32(onecon(:,:)),mapkeys),2);
%             typenum(fidx,type) = sum(typesel);
%             typesigcon=onecon(typesel,:);    
%             con_com_diff=arrayfun(@(x) com_str.(skey{1}).(samp)(x),int32(onecon(typesel,:)));
%             dirsel=(con_com_diff(:,2)-con_com_diff(:,1)>(50/250)) & (con_com_diff(:,2)-con_com_diff(:,1)<=8);
%             dirsigcon=typesigcon(dirsel,:);
%             sumcom(fidx,type) = size(dirsigcon,1);
%         end
%     end
% end
% 
% 
% 
% a = (sumcom(:,1) + sumcom(:,2))./(typenum(:,1)+typenum(:,2));
% b = (sumcom(:,3) + sumcom(:,4))./(typenum(:,3)+typenum(:,4));
% %pairs = chainsumCE(:,:,1);
% sel = all(typenum>=5,2) & ~isnan(a) & ~isnan(b);
% 
% meanCOMconnC = mean(a(sel));
% meanCOMconnE = mean(b(sel));
% 
% semCOMconnC = std(a(sel))/sqrt(sum(sel));
% semCOMconnE = std(b(sel))/sqrt(sum(sel));
% 
% fh=figure('Color','w','Position',[100,100,235,270]);
% hold on
% bar([meanCOMconnC,meanCOMconnE]*100,'FaceColor','w')
% errorbar(1:2,[meanCOMconnC,meanCOMconnE]*100,[semCOMconnC,semCOMconnE]*100,'k.','CapSize',12)
% set(gca,'XTick',1:2,'XTickLabel',{'Correct Trial','Error Trials'},'XTickLabelRotation',30)
% ylabel('% 2 Neuron Wave')