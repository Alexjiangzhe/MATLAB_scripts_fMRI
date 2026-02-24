function [] = singletrial_batch()

% Author: Maureen Ritchey, 10-2012
% Modified to auto-detect subjects and time folders
% Modified to detect existing outputs and only process missing data
% Modified to support manual subject and time window input

%% USER INFORMATION

% Directory for the existing model (主路径)
modeldir = '/Users/alexjiangzhe/Documents/My_project/1level_traning/output';

% Directory for the new model. If it does not exist, it will be created. (输出路径)
stdir = '/Users/alexjiangzhe/Documents/My_project/1level_traning/output/single_trial_model\';

% Flag for estimating model (0=no estimation, 1=estimation)
estimate = 1;

% If desired, specify a condition name to lump together (no indiv trial
% regressors). This is useful if the original SPM.mat model included a
% 'catch-all' condition of no interest. Set to 'NONE' if not desired.
lump_conditions = 'NONE'; % 'OTHERS'

% Type of model to run (1=multi-regressor approach, 2=multi-model approach)
modeltype = 1;

% Option to discard extra files in multi-model approach to save space
% (1=discard, 0=leave them). Always keeps regs and covs files.
discard_mm_files = 1;

% Flag for skipping already processed folders (1=skip, 0=reprocess/overwrite)
skip_existing = 1;

%% ========== 手动输入选项 ==========

% 被试选择模式：
%   'AUTO'     - 自动检测所有被试
%   'MANUAL'   - 手动指定被试
manual_mode = 'MANUAL';  % 改为 'MANUAL' 启用手动输入

% 手动指定被试（当 manual_mode = 'MANUAL' 时生效）
% 支持多种格式：
%   - 单个被试：{'1'}
%   - 多个被试：{'1', '2', '3', '10'}
%   - 范围：可以用数字数组，脚本会自动转换
manual_subjects = {'1'};  % 示例：被试1、2、3

% 手动指定时间窗（可选，留空则处理该被试的所有时间窗）
% 支持多种格式：
%   - 空数组：{}  - 处理所有时间窗
%   - 指定时间窗：{'time1', 'time2'}
manual_timewindows = {'0-1.5', '0-3', '1.5-3'};  % 示例：留空表示处理所有时间窗

% 注意：如果同时指定了被试和时间窗，将只处理指定被试的指定时间窗
% 如果只指定被试不指定时间窗，将处理指定被试的所有时间窗

%% ====================================

%% MAIN CODE

fprintf('\n╔════════════════════════════════════════════════════════════╗\n');
fprintf('║     Single Trial Batch Processing - Version 2.2          ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n');

% 决定使用自动检测还是手动输入
if strcmpi(manual_mode, 'MANUAL')
    fprintf('\n[手动模式] 使用用户指定的被试和时间窗\n');
    
    % 验证手动输入的被试
    input_subjects = validate_manual_subjects(manual_subjects, modeldir);
    
    if isempty(input_subjects)
        error('指定的被试在输入目录中不存在！请检查 manual_subjects 设置。');
    end
    
    fprintf('手动指定被试: %s\n', strjoin(input_subjects, ', '));
    
    if ~isempty(manual_timewindows)
        fprintf('手动指定时间窗: %s\n', strjoin(manual_timewindows, ', '));
    else
        fprintf('时间窗: 处理所有可用时间窗\n');
    end
else
    fprintf('\n[自动模式] 扫描输入目录中的所有被试\n');
    
    % 自动检测被试号
    fprintf('Input directory: %s\n', modeldir);
    subject_dirs = dir(modeldir);
    subject_dirs = subject_dirs([subject_dirs.isdir] & ~ismember({subject_dirs.name}, {'.', '..'}));
    
    % 提取数字被试号
    input_subjects = {};
    for i = 1:length(subject_dirs)
        if ~isempty(regexp(subject_dirs(i).name, '^\d+$', 'once'))  % 匹配纯数字文件夹名
            input_subjects{end+1} = subject_dirs(i).name;
        end
    end
    
    if isempty(input_subjects)
        error('No numeric subject folders found in %s', modeldir);
    end
    
    % 按数字顺序排序
    input_subjects = sort_numeric_subjects(input_subjects);
    fprintf('Found %d subjects: %s\n', length(input_subjects), strjoin(input_subjects, ', '));
end

% 检测输出路径中已有的被试和时间段
fprintf('\n=== STEP 1: Scanning output directory for existing results ===\n');
fprintf('Output directory: %s\n', stdir);
existing_data = struct();

if exist(stdir, 'dir')
    output_subject_dirs = dir(stdir);
    output_subject_dirs = output_subject_dirs([output_subject_dirs.isdir] & ~ismember({output_subject_dirs.name}, {'.', '..'}));
    
    for i = 1:length(output_subject_dirs)
        subj = output_subject_dirs(i).name;
        if ~isempty(regexp(subj, '^\d+$', 'once'))
            % 检测该被试下已有的时间段
            subj_output_path = fullfile(stdir, subj);
            time_dirs = dir(subj_output_path);
            time_dirs = time_dirs([time_dirs.isdir] & ~ismember({time_dirs.name}, {'.', '..'}));
            
            existing_times = {};
            for j = 1:length(time_dirs)
                time_folder = time_dirs(j).name;
                % 检查该时间段是否有完整的输出文件
                is_complete = false;
                
                if modeltype == 1
                    spm_file = fullfile(subj_output_path, time_folder, 'SPM.mat');
                    beta_file = fullfile(subj_output_path, time_folder, 'beta_info.mat');
                    if exist(spm_file, 'file') && exist(beta_file, 'file')
                        is_complete = true;
                    end
                elseif modeltype == 2
                    betadir = fullfile(subj_output_path, time_folder, 'betas');
                    beta_file = fullfile(betadir, [subj '_beta_info.mat']);
                    if exist(betadir, 'dir') && exist(beta_file, 'file')
                        beta_files = dir(fullfile(betadir, '*.img'));
                        if ~isempty(beta_files)
                            is_complete = true;
                        end
                    end
                end
                
                if is_complete
                    existing_times{end+1} = time_folder;
                end
            end
            
            if ~isempty(existing_times)
                % MATLAB字段名不能以数字开头，添加前缀
                field_name = ['sub_' subj];
                existing_data.(field_name) = existing_times;
                fprintf('  Subject %s: found %d complete time folders: %s\n', ...
                    subj, length(existing_times), strjoin(existing_times, ', '));
            end
        end
    end
    
    if isempty(fieldnames(existing_data))
        fprintf('No existing complete results found in output directory.\n');
    end
else
    fprintf('Output directory does not exist yet. Will process all subjects.\n');
end

% 统计信息
fprintf('\n=== STEP 2: Planning processing tasks ===\n');
total_to_process = 0;
total_to_skip = 0;
processing_plan = struct();

for iSub = 1:length(input_subjects)
    subj = input_subjects{iSub};
    
    % 获取输入目录中该被试的所有时间段
    subject_input_path = fullfile(modeldir, subj);
    time_dirs = dir(subject_input_path);
    time_dirs = time_dirs([time_dirs.isdir] & ~ismember({time_dirs.name}, {'.', '..'}));
    
    if isempty(time_dirs)
        fprintf('Warning: No time folders found for subject %s\n', subj);
        continue;
    end
    
    all_input_times = {time_dirs.name};
    
    % 如果手动模式且指定了时间窗，则筛选
    if strcmpi(manual_mode, 'MANUAL') && ~isempty(manual_timewindows)
        % 只保留手动指定的时间窗
        input_times = intersect(all_input_times, manual_timewindows);
        
        if isempty(input_times)
            fprintf('Warning: Subject %s has no specified time windows (%s)\n', ...
                subj, strjoin(manual_timewindows, ', '));
            continue;
        end
    else
        input_times = all_input_times;
    end
    
    % 确定需要处理的时间段
    % MATLAB字段名不能以数字开头，使用前缀
    field_name = ['sub_' subj];
    if skip_existing && isfield(existing_data, field_name)
        existing_times = existing_data.(field_name);
        times_to_process = setdiff(input_times, existing_times);
        times_to_skip = intersect(input_times, existing_times);
    else
        times_to_process = input_times;
        times_to_skip = {};
    end
    
    % 保存处理计划（使用带前缀的字段名）
    processing_plan.(field_name).to_process = times_to_process;
    processing_plan.(field_name).to_skip = times_to_skip;
    
    total_to_process = total_to_process + length(times_to_process);
    total_to_skip = total_to_skip + length(times_to_skip);
    
    if ~isempty(times_to_process)
        fprintf('Subject %s: will process %d time folders: %s\n', ...
            subj, length(times_to_process), strjoin(times_to_process, ', '));
    end
    if ~isempty(times_to_skip)
        fprintf('Subject %s: will skip %d time folders: %s\n', ...
            subj, length(times_to_skip), strjoin(times_to_skip, ', '));
    end
end

fprintf('\n>>> Total time folders to process: %d\n', total_to_process);
fprintf('>>> Total time folders to skip: %d\n', total_to_skip);

if total_to_process == 0
    fprintf('\n>>> All specified data already processed! Nothing to do.\n');
    return;
end

% 开始处理
fprintf('\n=== STEP 3: Processing data ===\n');
failed_list = {};
successful_count = 0;

for iSub = 1:length(input_subjects)
    subj = input_subjects{iSub};
    field_name = ['sub_' subj];  % 添加前缀
    
    if ~isfield(processing_plan, field_name) || isempty(processing_plan.(field_name).to_process)
        continue;
    end
    
    fprintf('\n=== Processing Subject %s ===\n', subj);
    times_to_process = processing_plan.(field_name).to_process;
    
    for iTime = 1:length(times_to_process)
        time_folder = times_to_process{iTime};
        fprintf('\n--- Processing time folder: %s ---\n', time_folder);
        
        try
            % 加载SPM文件
            spm_file_path = fullfile(modeldir, subj, time_folder, 'SPM.mat');
            fprintf('Loading SPM file: %s\n', spm_file_path);
            
            if ~exist(spm_file_path, 'file')
                fprintf('Warning: SPM.mat not found in %s\n', spm_file_path);
                failed_list{end+1} = sprintf('%s/%s', subj, time_folder);
                continue;
            end
            
            load(spm_file_path);
            
            % 创建输出目录
            outputdir = fullfile(stdir, subj, time_folder, filesep);
            if ~exist(outputdir, 'dir')
                fprintf('Creating directory: %s\n', outputdir);
                mkdir(outputdir);
            end
            
            % 获取模型信息
            fprintf('Getting model information...\n');
            files = SPM.xY.P;
            fprintf('Modeling %i timepoints across %i sessions.\n', size(files,1), length(SPM.Sess));
            
            % MULTI-REGRESSOR APPROACH
            if modeltype == 1
                
                % 设置beta信息
                trialinfo = {'beta_number' 'session' 'condition' 'condition_rep' 'number_onsets' 'first_onset' 'beta_name'};
                counter = 1;
                
                % 循环处理sessions
                for iSess=1:length(SPM.Sess)
                    rows = SPM.Sess(iSess).row;
                    sess_files = files(rows',:);
                    sess_files = cellstr(sess_files);
                    covariates = SPM.Sess(iSess).C.C;
                    
                    onsets = {};
                    durations = {};
                    names = {};
                    
                    for j=1:length(SPM.Sess(iSess).U)
                        % 检查特殊条件名称
                        if strfind(SPM.Sess(iSess).U(j).name{1},lump_conditions)
                            onsets = [onsets SPM.Sess(iSess).U(j).ons'];
                            durations = [durations SPM.Sess(iSess).U(j).dur'];
                            cur_name = [SPM.Sess(iSess).U(j).name{1}];
                            names = [names cur_name];
                            curinfo = {counter iSess SPM.Sess(iSess).U(j).name{1} [1] length(SPM.Sess(iSess).U(j).ons) SPM.Sess(iSess).U(j).ons(1) cur_name};
                            trialinfo = [trialinfo; curinfo];
                            counter = counter + 1;
                        else
                            % 为每个个体试次设置回归器
                            for k=1:length(SPM.Sess(iSess).U(j).ons)
                                onsets = [onsets SPM.Sess(iSess).U(j).ons(k)];
                                durations = [durations SPM.Sess(iSess).U(j).dur(k)];
                                cur_name = [SPM.Sess(iSess).U(j).name{1} '_' num2str(k)];
                                names = [names cur_name];
                                curinfo = {counter iSess SPM.Sess(iSess).U(j).name{1} k length(SPM.Sess(iSess).U(j).ons(k)) SPM.Sess(iSess).U(j).ons(k) cur_name};
                                trialinfo = [trialinfo; curinfo];
                                counter = counter + 1;
                            end
                        end
                    end
                    
                    % 保存回归器onset文件
                    fprintf('Saving regressor onset files for Session %i: %i trials included\n',iSess,length(names));
                    regfile = [outputdir 'st_regs_run' num2str(iSess) '.mat'];
                    save(regfile,'names','onsets','durations');
                    
                    % 保存协变量
                    covfile = [outputdir 'st_covs_run' num2str(iSess) '.txt'];
                    dlmwrite(covfile,covariates,'\t');
                    if ~isempty(covariates)
                        for icov = 1:size(covariates,2)
                            curinfo = {counter iSess 'covariate' icov 1 0 strcat('covariate',num2str(icov))};
                            trialinfo = [trialinfo; curinfo];
                            counter = counter + 1;
                        end
                    end
                    
                    % 创建matlabbatch
                    if iSess==1
                        matlabbatch = create_spm_init(outputdir,SPM);
                    end
                    matlabbatch = create_spm_sess(matlabbatch,iSess,sess_files,regfile,covfile,SPM);
                    
                end
                
                % 运行matlabbatch创建SPM.mat文件
                fprintf('Creating SPM.mat file: %s\n',[outputdir 'SPM.mat']);
                spm_jobman('initcfg')
                spm('defaults', 'FMRI');
                spm_jobman('serial', matlabbatch);
                
                if estimate > 0
                    fprintf('Estimating model from SPM.mat file.\n');
                    spmfile = [outputdir 'SPM.mat'];
                    matlabbatch = estimate_spm(spmfile);
                    spm_jobman('serial', matlabbatch);
                end
                
                % 保存beta信息
                infofile = [outputdir 'beta_info.mat'];
                save(infofile,'trialinfo');
                
                clear SPM matlabbatch
                
            elseif modeltype == 2
                
                % MULTI-MODEL APPROACH
                betadir = [outputdir 'betas/'];
                if ~exist(betadir,'dir')
                    mkdir(betadir)
                end
                
                trialinfo = {'beta_number' 'session' 'condition' 'condition_rep' 'number_onsets' 'beta_name' 'trial_dir' 'beta_name'};
                counter = 1;
                
                for iSess=1:length(SPM.Sess)
                    rows = SPM.Sess(iSess).row;
                    sess_files = files(rows',:);
                    sess_files = cellstr(sess_files);
                    covariates = SPM.Sess(iSess).C.C;
                    
                    for j=1:length(SPM.Sess(iSess).U)
                        if strfind(SPM.Sess(iSess).U(j).name{1},lump_conditions)
                            if strcmp(lump_conditions,'NONE')
                            else
                                fprintf('Ignoring condition: %s\n',lump_conditions);
                            end
                        else
                            for k=1:length(SPM.Sess(iSess).U(j).ons)
                                
                                onsets = {};
                                durations = {};
                                names = {};
                                onsets = [onsets SPM.Sess(iSess).U(j).ons(k)];
                                durations = [durations SPM.Sess(iSess).U(j).dur(k)];
                                cur_name = [SPM.Sess(iSess).U(j).name{1} '_' num2str(k)];
                                names = [names cur_name];
                                
                                otheronsets = [];
                                otherdurations = [];
                                othernames = {};
                                for jj=1:length(SPM.Sess(iSess).U)
                                    for kk=1:length(SPM.Sess(iSess).U(jj).ons)
                                        if j~=jj & k~=kk
                                            otheronsets = [otheronsets SPM.Sess(iSess).U(jj).ons(kk)];
                                            otherdurations = [otherdurations SPM.Sess(iSess).U(jj).dur(kk)];
                                            othernames = ['OTHERTRIALS'];
                                        end
                                    end
                                end
                                
                                onsets = [onsets otheronsets];
                                durations = [durations otherdurations];
                                names = [names othernames];
                                
                                trialdir = [outputdir 'Sess' num2str(iSess) '/' cur_name '/'];
                                if ~exist(trialdir,'dir')
                                    mkdir(trialdir)
                                end
                                
                                curinfo = {counter iSess SPM.Sess(iSess).U(j).name{1} k length(SPM.Sess(iSess).U(j).ons(k)) cur_name trialdir ['Sess' num2str(iSess) '_' cur_name '.img']};
                                trialinfo = [trialinfo; curinfo];
                                counter = counter + 1;
                                
                                regfile = [trialdir 'st_regs.mat'];
                                save(regfile,'names','onsets','durations');
                                
                                covfile = [trialdir 'st_covs.txt'];
                                dlmwrite(covfile,covariates,'\t');
                                
                                matlabbatch = create_spm_init(trialdir,SPM);
                                matlabbatch = create_spm_sess(matlabbatch,1,sess_files,regfile,covfile,SPM);
                                
                                fprintf('Creating SPM.mat file: %s\n',[trialdir 'SPM.mat']);
                                if counter == 1
                                    spm_jobman('initcfg')
                                    spm('defaults', 'FMRI');
                                end
                                spm_jobman('serial', matlabbatch);
                                
                                if estimate
                                    fprintf('Estimating model from SPM.mat file.\n');
                                    spmfile = [trialdir 'SPM.mat'];
                                    matlabbatch = estimate_spm(spmfile);
                                    spm_jobman('serial', matlabbatch);
                                    clear matlabbatch
                                    
                                    copyfile([trialdir 'beta_0001.img'],[betadir 'Sess' num2str(iSess) '_' cur_name '.img']);
                                    copyfile([trialdir 'beta_0001.hdr'],[betadir 'Sess' num2str(iSess) '_' cur_name '.hdr']);
                                    
                                    if discard_mm_files
                                        prevdir = pwd;
                                        cd(trialdir);
                                        delete SPM*; delete *.hdr; delete *.img;
                                        cd(prevdir);
                                    end
                                    
                                end
                            end
                            
                        end
                    end
                end
                
                infofile = [betadir subj '_beta_info.mat'];
                save(infofile,'trialinfo');
                
            else
                error('Specify model type as 1 or 2');
            end
            
            clear SPM
            fprintf('✓ Successfully processed %s - %s\n', subj, time_folder);
            successful_count = successful_count + 1;
            
        catch ME
            fprintf('✗ Error processing %s - %s: %s\n', subj, time_folder, ME.message);
            failed_list{end+1} = sprintf('%s/%s', subj, time_folder);
            continue;
        end
    end
end

% 最终总结
fprintf('\n========================================\n');
fprintf('       PROCESSING SUMMARY\n');
fprintf('========================================\n');
if strcmpi(manual_mode, 'MANUAL')
    fprintf('Mode: MANUAL\n');
    fprintf('Specified subjects: %d\n', length(input_subjects));
else
    fprintf('Mode: AUTO\n');
    fprintf('Total subjects scanned: %d\n', length(input_subjects));
end
fprintf('Total time folders in scope: %d\n', total_to_process + total_to_skip);
fprintf('Skipped (already complete): %d\n', total_to_skip);
fprintf('Attempted to process: %d\n', total_to_process);
fprintf('Successfully processed: %d\n', successful_count);
fprintf('Failed: %d\n', length(failed_list));

if ~isempty(failed_list)
    fprintf('\nFailed items:\n');
    for i = 1:length(failed_list)
        fprintf('  - %s\n', failed_list{i});
    end
end

fprintf('========================================\n');

end

%% HELPER FUNCTIONS

function valid_subjects = validate_manual_subjects(manual_subjects, modeldir)
% 验证手动输入的被试是否存在于输入目录中
valid_subjects = {};

for i = 1:length(manual_subjects)
    subj = manual_subjects{i};
    subj_path = fullfile(modeldir, subj);
    
    if exist(subj_path, 'dir')
        valid_subjects{end+1} = subj;
    else
        fprintf('Warning: Subject %s not found in input directory, skipping.\n', subj);
    end
end
end

function sorted_subjects = sort_numeric_subjects(subjects)
% 将字符串数字转换为数值进行排序
numeric_values = cellfun(@str2double, subjects);
[~, idx] = sort(numeric_values);
sorted_subjects = subjects(idx);
end

%% SUBFUNCTIONS

function [matlabbatch] = create_spm_init(outputdir,SPM)
% subfunction for initializing the matlabbatch structure to create the SPM
matlabbatch{1}.spm.stats.fmri_spec.dir = {outputdir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = SPM.xBF.UNITS;
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = SPM.xY.RT;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = SPM.xBF.T;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = SPM.xBF.T0;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = SPM.xBF.Volterra;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
if isempty(SPM.xM.VM)
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
else
    matlabbatch{1}.spm.stats.fmri_spec.mask = {SPM.xM.VM.fname};
end
matlabbatch{1}.spm.stats.fmri_spec.cvi = SPM.xVi.form;

end

function [matlabbatch] = create_spm_sess(matlabbatch,iSess,sess_files,regfile,covfile,SPM)
% subfunction for adding sessions to the matlabbatch structure
matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).scans = sess_files;
matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).multi = {regfile};
matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).multi_reg = {covfile};
matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).hpf = SPM.xX.K(iSess).HParam;

end

function [matlabbatch] = estimate_spm(spmfile)
% subfunction for creating a matlabbatch structure to estimate the SPM
matlabbatch{1}.spm.stats.fmri_est.spmmat = {spmfile};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
end