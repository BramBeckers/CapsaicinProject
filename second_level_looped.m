% Second-level analysis (pseudo-block pharmacological MRI)
%Authors: Lukas van Oudenhove, Jie Wu, KU Leuven

%root = ''; specify folder with first level results
cd(root);spm('defaults','FMRI');%set default
%file = dir(''); specify name of first level results folders

%mkdir(''); create second level folder
sbjn = 18;      %subjectnr
binn = 45;      %amount of timebins
%matlabbatch{1}.spm.stats.factorial_design.dir = {''}; input created folder

for i=1:length(file);
cd([root file(i).name]);
file_con = spm_select('FPList', [root file(i).name], ['^','con_*.*\.nii']);   %character array
filecon{i} = cellstr(file_con);
j = filecon{i} (1:45,1);                                                      
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(i).scans = j;

matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(i).conds = [1:binn];
cd(root);
end;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'Interaction effect';
matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = {diff(eye(binn))};
matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.fcon.name = 'secondary';
matlabbatch{3}.spm.stats.con.consess{2}.fcon.weights = {eye(binn) ones(binn,sbjn)/sbjn};
matlabbatch{3}.spm.stats.con.consess{2}.fcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
cd '';
eval('save Batch_2nd_level_total.mat matlabbatch');   %save batch job and rename
cd(root);
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);%run batch
clear matlabbatch;%clear the batch for next
spm('defaults','FMRI');%set default