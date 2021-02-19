% First-level analysis (pseudo-block pharmacological MRI)
%Authors: Lukas van Oudenhove, Jie Wu, KU Leuven

% adapted by Bram Beckers, Maastricht University
% date:      May 2019
%_________________________________________________________________________
close all hidden
clear all hidden

% define settings of SPMl
spm('defaults','FMRI')

% define directory where denoised images are (one single directory, move/copy
% images there prior to running the script if needed)
% datadir = ''; specify path to input folder
cd(datadir);

% sub_fructans = dir('niftiDATA_subjcon*_fructans.nii');
% sub_capsaicin = dir(''); specify filenames of condition A, eg
% data_condA_subj*.nii (use asterisk to loop through multiple files)
%sub_saline = dir(''); specify filename of condition B, eg
% data_condB_subj*.nii
% template_struct_image = 'wc1.*\.nii$';
%-------------------------------------------------------------
% nr_subjects = size(sub_capsaicin,1);
% loop creation of 2 batches
for i = 1:length(sub_capsaicin) %% not really sure if this name is correct
    clear matlabbatch
    matlabbatch = struct([]); % create empty structures
    analysis = ['first_level_',sub_capsaicin(i).name(11:20)]; %% what does the 11:20 specify (related to analysis params?)?
    mkdir(analysis); %create new analysis folder, renamed!
%     file_fructans = spm_select('Extlist', datadir, sub_fructans(i).name,1:1416);
    file_capsaicin = spm_select('Extlist', datadir, sub_capsaicin(i).name,1:1350);
    file_saline = spm_select('Extlist', datadir, sub_saline(i).name,1:1350);
    cd(analysis);

    % FIRST BATCH: MODEL SPECIFICATION
    %--------------------------------
    % make matlab batch structure
    % ---------------------------
    matlabbatch{1}.spm.stats.fmri_spec.dir = {analysis};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.6;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 50;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 25; %%Despite multiband acquisition microtime onset can be kept at middle slice since slice time correction has been applied
    % SESSION 1 Capsaicin
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(strcat(datadir, file_capsaicin));
    for z = 1:45  %loop create timebins for first session
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(z).name = ['capsaicin_bin' num2str(z)];
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(z).onset = 25 * z + 201; %%217 + 24 is 241 (number of baseline volumes)
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(z).duration = 25;
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(z).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(z).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(z).orth = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 2160; %%(Could it be that this number has to be increased in our case? 'infinite'? --> volumes*TR)
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
    end;
    
    % SESSION 2 saline
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = cellstr(strcat(datadir, file_saline));
    for y = 1:45  %loop create timebins for second session
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(y).name = ['saline_bin' num2str(y)];
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(y).onset = 25 * y + 201;  
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(y).duration = 25;
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(y).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(y).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(y).orth = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 2160;
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = {''};
    end;

%     % SESSION 3 (extra, not needed, just in case you have more conditions)
%     matlabbatch{1}.spm.stats.fmri_spec.sess(3).scans = cellstr(strcat(datadir, file_saline));
%     for x = 1:45  %loop create timebins for third session
%     matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(x).name = ['saline_bin' num2str(x)];
%     matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(x).onset = 25 * x + 201;  
%     matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(x).duration = 25;
%     matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(x).tmod = 0;
%     matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(x).pmod = struct('name', {}, 'param', {}, 'poly', {});
%     matlabbatch{1}.spm.stats.fmri_spec.sess(3).cond(x).orth = 1;
%     matlabbatch{1}.spm.stats.fmri_spec.sess(3).hpf = 3540;
%     matlabbatch{1}.spm.stats.fmri_spec.sess(3).multi = {''};
%     end;
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8; 

    % SECOND BATCH: ESTIMATION
    %------------------------
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    % THIRD BATCH: CONTRASTS
    %----------------------
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    for a = 1:45 %loop for contrast manager 
    matlabbatch{3}.spm.stats.con.consess{a}.tcon.name = ['capsaicin_saline_bin' num2str(a)];
    matlabbatch{3}.spm.stats.con.consess{a}.tcon.weights = [zeros(1, a-1) 1 zeros(1,44) -1];   %first 1 indicates amount of rows ( but only 1 row here, since it is repeated for each timebin seperate ), a-1 indicates amount of colums. 
    matlabbatch{3}.spm.stats.con.consess{a}.tcon.sessrep = 'none';
    end;
%     for b = 1:49 %loop for contrast manager
%         d = b + 49;
%     matlabbatch{3}.spm.stats.con.consess{d}.tcon.name = ['fructans_saline_bin' num2str(b)];
%     matlabbatch{3}.spm.stats.con.consess{d}.tcon.weights = [zeros(1, b-1) 1 zeros(1,97) -1];   %first 1 indicates amount of rows ( but only 1 row here, since it is repeated for each timebin seperate ), b-1 indicates amount of colums. 
%     matlabbatch{3}.spm.stats.con.consess{d}.tcon.sessrep = 'none';
%     end;
%      for c = 1:49 %loop for contrast manager
%          e = c + 98;
%     matlabbatch{3}.spm.stats.con.consess{e}.tcon.name = ['glucose_saline_bin' num2str(c)];
%     matlabbatch{3}.spm.stats.con.consess{e}.tcon.weights = [zeros(1, 49) zeros(1,c-1) 1 zeros(1,48) -1];   %first 1 indicates amount of rows ( but only 1 row here, since it is repeated for each timebin seperate ), c-1 indicates amount of colums. 
%     matlabbatch{3}.spm.stats.con.consess{e}.tcon.sessrep = 'none';
%     end;
   matlabbatch{3}.spm.stats.con.delete = 0;
    % SAVE BATCHES AND RUN
    %---------------------
    eval(['save design_first_level_' sub_capsaicin(i).name(11:20) '.mat matlabbatch']); 
    cd(datadir)
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
end;