function countWaves(opt)
% prepare functional coupling linked neuronal-wave data
arguments
    opt.peak (1,1) logical = false
    opt.strict (1,1) logical = true %strict FC criteria
%     opt.min_diff (1,1)
    opt.screen (1,1) logical = false
end
load('sums_conn.mat','sums_conn_str');
%chainsum = cell(length(sums_conn_str),1);
for fidx=1:length(sums_conn_str)
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
    onecom=wave.get_com_map('onepath',strrep(sums_conn_str(fidx).folder,'zx/neupix','xd/data'),'peak',opt.peak,'pathid',fidx); %per su center of mass, normalized FR
end
        
        
        