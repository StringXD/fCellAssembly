function com_str_=get_com_map(opt)
arguments
    opt.onepath (1,:) char = ''
    opt.peak (1,1) logical = false
    opt.curve (1,1) logical = false
    opt.per_sec_stats (1,1) logical = false
    opt.selidx_curve (1,1) logical = false
    opt.pathid (1,1)
end
persistent com_str onepath_
if isempty(onepath_), onepath_='';end
if isempty(com_str) || ~strcmp(opt.onepath, onepath_)
    meta_str=ephys.util.load_meta('type','neupix');
    %homedir=ephys.util.getHomedir('type','raw');
    com_str=struct();
    
    
        %dpath=regexp(opt.onepath,'(?<=SPKINFO[\\/]).*$','match','once');
        %fpath=fullfile(homedir,dpath,'FR_All_1000.hdf5');
        dpath = opt.onepath;
        fpath=fullfile(dpath,'FR_All_250.hdf5');
        %pc_stem=replace(dpath,'/',filesep());
        t=strsplit(dpath,'/');
        pc_stem=[t{end-1},'\',t{end}];
        sesssel=startsWith(meta_str.allpath,pc_stem);
        %save('sesssel.mat','meta_str','pc_stem','-v7.3');
        
        fr=h5read(fpath,'/FR_All');
        trial=h5read(fpath,'/Trials');
        suid=h5read(fpath,'/SU_id');
        nmcid=meta_str.allcid(meta_str.mem_type==0 & sesssel.');
        nm=find(ismember(suid,nmcid));
        mcid1=meta_str.allcid(meta_str.mem_type==2 & sesssel.');
        msel1=find(ismember(suid,mcid1));
        mcid2=meta_str.allcid(meta_str.mem_type==4 & sesssel.');
        msel2=find(ismember(suid,mcid2));
        if ~isempty(msel1) && ~isempty(msel2)
        %opt.pathid=ephys.path2opt.pathid(pc_stem);
        % suppose non-mem neuron com follow uniform distribution
        com_str.(['s',num2str(opt.pathid)]).s0=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(opt.pathid)]).s1=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(opt.pathid)]).s2=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(opt.pathid)]).s1e=containers.Map('KeyType','int32','ValueType','any'); % error trials
        com_str.(['s',num2str(opt.pathid)]).s2e=containers.Map('KeyType','int32','ValueType','any'); % error trials
        com_str.(['s',num2str(opt.pathid)]).s1curve=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(opt.pathid)]).s2curve=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(opt.pathid)]).s1ecurve=containers.Map('KeyType','int32','ValueType','any');
        com_str.(['s',num2str(opt.pathid)]).s2ecurve=containers.Map('KeyType','int32','ValueType','any');
        %     if sum(trial(:,9))<40,continue;end %meta data obtained from processed
        %     welltrained dataset
        s1sel=trial(:,5)==4 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;
        s2sel=trial(:,5)==8 & trial(:,8)==6 & trial(:,9)>0 & trial(:,10)>0;
        s1err=trial(:,5)==4 & trial(:,8)==6 & trial(:,10)==0;
        s2err=trial(:,5)==8 & trial(:,8)==6 & trial(:,10)==0;
        sess=['s',num2str(opt.pathid)];
        com_str=per_su_process(sess,suid,nm,fr,s1sel|s2sel,0,com_str,'s0',opt);
        com_str=per_su_process(sess,suid,msel1,fr,s1sel,s2sel,com_str,'s1',opt);
        com_str=per_su_process(sess,suid,msel2,fr,s2sel,s1sel,com_str,'s2',opt);
        if sum(s1err)>=5 && sum(s2err)>=5
            com_str=per_su_process(sess,suid,msel1,fr,s1err,s2err,com_str,'s1e',opt);
            com_str=per_su_process(sess,suid,msel2,fr,s2err,s1err,com_str,'s2e',opt);
        end
        end
end

com_str_=com_str;
save(sprintf('com_str_Newcrit_%d.mat', opt.pathid),'com_str');
end

function com_str=per_su_process(sess,suid,msel,fr,pref_sel,nonpref_sel,com_str,samp,opt)
if strcmp(samp,'s0')
    for su=reshape(msel,1,[])
        basemm=mean(mean(mean(squeeze(fr(pref_sel,su,17:40)))));
        mm=smooth(squeeze(mean(fr(pref_sel,su,:))),5).';
        mm_pref=mm(17:40)-basemm;
        mm_pref=mm_pref./max(mm_pref);
        mm_pref(mm_pref<0)=0;
        if opt.peak
            [~,pidx]=max(mm_pref);
            com_str.(sess).(samp)(suid(su))=pidx;
        else
            com=sum((1:24).*mm_pref)./sum(mm_pref);
            com_str.(sess).(samp)(suid(su))=com;
        end
    end

else
    for su=reshape(msel,1,[])
        basemm=mean([mean(squeeze(fr(pref_sel,su,17:40)));mean(squeeze(fr(nonpref_sel,su,17:40)))]);
        if ~opt.per_sec_stats
            basemm=mean(basemm);
        end
        mm=smooth(squeeze(mean(fr(pref_sel,su,:))),5).';
        mm_pref=mm(17:40)-basemm;
        mm_pref=mm_pref./max(mm_pref);
        mm_pref(mm_pref<0)=0;
        if opt.peak
            [~,pidx]=max(mm_pref);
            com_str.(sess).(samp)(suid(su))=pidx;
        else
            com=sum((1:24).*mm_pref)./sum(mm_pref);
            com_str.(sess).(samp)(suid(su))=com;
        end
%     if opt.selidx_curve
%         fr_pref=squeeze(mean(fr(pref_sel,su,17:40)));
%         fr_nonp=squeeze(mean(fr(nonpref_sel,su,17:40)));
%         com_str.(sess).([samp,'curve'])(suid(su))...
%             =(fr_pref-fr_nonp)./(fr_pref+fr_nonp);
%     else
%         com_str.(sess).([samp,'curve'])(suid(su))...
%             =mm_pref;
%     end
    end
end

end