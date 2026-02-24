clear; clc;
spm('defaults', 'FMRI');
spm_jobman('initcfg');

%% ===== 用户配置区域（按需要改路径） =====
% 需要跑的被试
subjects_to_run = {'1'};
%'25','26','27','23','24','28','29','30',
% 需要建模的持续时间窗口（对应不同的列）
durations_to_run = {'0-6','1.5-3'};
%,'0-6','3-6'
% 条件名称，要和每个被试 ONSET 文件夹里的 xlsx 文件名完全对应
% 例如 D:\TIAOJIANJIQIXUEXI\ONSET\16\distancing.xlsx
conditions = {'blank','distancing','distraction','inhibition','nature'};

% 路径
onset_base_dir   = '/Users/alexjiangzhe/Documents/My_project/1level_traning/ONSET';    % 每个被试一个文件夹，里面是5个xlsx
preproc_base_dir = '//Users/alexjiangzhe/Documents/My_project/1level_traning';      % 例如 D:\tiaojianjiqixuexi\pre
output_base_dir  = '/Users/alexjiangzhe/Documents/My_project/1level_traning/output'; % 结果输出路径（你可以改成自己想要的）

%% ===== 主循环 =====
for i = 1:length(subjects_to_run)
    subject_number = subjects_to_run{i};
    fprintf('====== 开始处理被试: %s ======\n', subject_number);

    try
        for j = 1:length(durations_to_run)
            duration_label = durations_to_run{j};
            fprintf('--- 分析时长: %s ---\n', duration_label);

            % ==== 根据持续时间选择 onset 所在列 和 duration ====
            % 当前 xlsx 列结构：
            % A: prepare.OffsetTim
            % B: Slide2.OnsetTim
            % C: slide2_0
            % D: slide2_1.5
            % E: slide2_3
            % F: slide2_4.5
            switch duration_label
                case '0-1.5'
                    % 从 0s 开始，持续 1.5s
                    onset_col_idx = 3;   % slide2_0
                    event_duration = 1.5;
                case '0-3'
                    % 从 0s 开始，持续 3s
                    onset_col_idx = 3;   % slide2_0
                    event_duration = 3;
                case '0-6'
                    % 从 0s 开始，持续 6s
                    onset_col_idx = 3;   % slide2_0
                    event_duration = 6;
                case '1.5-3'
                    % 从 1.5s 开始，持续 1.5s
                    onset_col_idx = 4;   % slide2_1.5
                    event_duration = 1.5;
                case '3-4.5'
                    % 从 3s 开始，持续 1.5s
                    onset_col_idx = 5;   % slide2_3
                    event_duration = 1.5;
                case '3-6'
                    % 从 3s 开始，持续 3s
                    onset_col_idx = 5;   % slide2_3
                    event_duration = 3;
                case '4.5-6'
                    % 从 4.5s 开始，持续 1.5s
                    onset_col_idx = 6;   % slide2_4.5
                    event_duration = 1.5;
                otherwise
                    error('未知时长标签: %s', duration_label);
            end



            % 创建输出目录：firstlevel\16\0-1.5\ 这样的结构
            output_dir = fullfile(output_base_dir, subject_number, duration_label);
            if ~exist(output_dir, 'dir'), mkdir(output_dir); end

            %% ===== 初始化 batch =====
            clear matlabbatch
            matlabbatch{1}.spm.stats.fmri_spec.dir = {output_dir};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.5;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 72;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 36;

            % -------------------- 只有 1 个 session --------------------
            sess = 1;

            % 被试预处理目录，例如 D:\tiaojianjiqixuexi\pre\16
            subj_pre_dir = fullfile(preproc_base_dir, subject_number);
            swuaf_dir    = fullfile(subj_pre_dir, 'swuaf');
            if ~exist(swuaf_dir, 'dir')
                error('未找到 swuaf 目录: %s', swuaf_dir);
            end

            % 搜索 3D NIfTI（不再硬性要求 1350 张，只要找到就行）
            nii_struct = dir(fullfile(swuaf_dir, '*.nii'));
            if isempty(nii_struct)
                error('被试 %s 未在 %s 找到任何 NIfTI 文件', subject_number, swuaf_dir);
            end

            % 按文件名排序并组装 scans
            [~, order] = sort({nii_struct.name});
            nii_struct = nii_struct(order);
            scans = cell(numel(nii_struct), 1);
            for ii = 1:numel(nii_struct)
                scans{ii,1} = fullfile(swuaf_dir, nii_struct(ii).name);
            end
            fprintf('    找到 %d 个 volume。\n', numel(nii_struct));
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess).scans = scans;

            % -------------------- 为每个 condition 加载 onset --------------------
            onset_dir = fullfile(onset_base_dir, subject_number);
            if ~exist(onset_dir, 'dir')
                error('未找到该被试的 ONSET 目录: %s', onset_dir);
            end

            for c = 1:length(conditions)
                cond_name = conditions{c};

                % 对应的 xlsx 文件，例如 D:\...\ONSET\16\distancing.xlsx
                onset_file = fullfile(onset_dir, [cond_name '.xlsx']);
                if ~exist(onset_file, 'file')
                    error('被试 %s 缺少 onset 文件: %s', subject_number, onset_file);
                end

                % 读取表（默认第一张 Sheet）
                T = readtable(onset_file);

                % 用列索引抽取对应列的数据（第 3~6 列）
                if width(T) < onset_col_idx
                    error('文件 %s 中列数不足，期望至少 %d 列', onset_file, onset_col_idx);
                end
                onset_data = T{:, onset_col_idx};
                onset_data = onset_data(~isnan(onset_data));
                onset_data = onset_data(:)';  % 行向量

                fprintf('    [✓] %s: %s 读取到 %d 个 onset（列索引=%d，dur=%.1f）\n', ...
                        subject_number, cond_name, numel(onset_data), onset_col_idx, event_duration);

                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).name     = cond_name;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).onset    = onset_data;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).duration = event_duration;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).tmod     = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess).cond(c).orth     = 1;
            end

            % 其它 session 级设置
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess).multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(sess).hpf   = 128;

            % ---- 运动参数：优先找 art_*.mat，其次 rp_*.txt ----
            rp_mat = dir(fullfile(subj_pre_dir, 'art_*.mat'));
            if isempty(rp_mat)
                rp_mat = dir(fullfile(swuaf_dir, 'art_*.mat'));
            end
            rp_txt = dir(fullfile(subj_pre_dir, 'rp_*.txt'));
            if isempty(rp_txt)
                rp_txt = dir(fullfile(swuaf_dir, 'rp_*.txt'));
            end

            if ~isempty(rp_mat)
                motion_file = fullfile(rp_mat(1).folder, rp_mat(1).name);
                fprintf('    [√] 使用MAT运动文件: %s\n', motion_file);
            elseif ~isempty(rp_txt)
                motion_file = fullfile(rp_txt(1).folder, rp_txt(1).name);
                fprintf('    [√] 使用TXT运动文件: %s\n', motion_file);
            else
                error('未找到运动文件 (art_* 或 rp_*): %s', subj_pre_dir);
            end

            matlabbatch{1}.spm.stats.fmri_spec.sess(sess).multi_reg = {motion_file};

            %% ===== 模型设置 =====
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt   = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
            matlabbatch{1}.spm.stats.fmri_spec.mask   = {''};
            matlabbatch{1}.spm.stats.fmri_spec.cvi    = 'AR(1)';

            %% ===== 模型估计 =====
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep( ...
                'fMRI model specification: SPM.mat File', ...
                substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

            %% ===== 对比设置（5 个条件各自一个 t-contrast） =====
            matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep( ...
                'Model estimation: SPM.mat File', ...
                substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                substruct('.','spmmat'));

            % 5 条主效应对比
            cons = {
                'blank',       [1 0 0 0 0]
                'distancing',  [0 1 0 0 0]
                'distraction', [0 0 1 0 0]
                'inhibition',  [0 0 0 1 0]
                'nature',      [0 0 0 0 1]
            };

            for k = 1:size(cons,1)
                matlabbatch{3}.spm.stats.con.consess{k}.tcon.name    = cons{k,1};
                matlabbatch{3}.spm.stats.con.consess{k}.tcon.weights = cons{k,2};
                matlabbatch{3}.spm.stats.con.consess{k}.tcon.sessrep = 'replsc';
            end
            matlabbatch{3}.spm.stats.con.delete = 0;

            %% ===== 保存并运行 batch =====
            batch_dir = fullfile(output_dir, 'batch');
            if ~exist(batch_dir, 'dir'), mkdir(batch_dir); end
            batch_filename = fullfile(batch_dir, sprintf('1st_batch_sub-%s_dur-%s.mat', subject_number, duration_label));
            save(batch_filename, 'matlabbatch');
            fprintf('  [已保存] 批处理文件至 %s\n', batch_filename);

            fprintf('>> 正在运行SPM Job: sub-%s, dur-%s\n', subject_number, duration_label);
            spm_jobman('run', matlabbatch);
            fprintf('>> 完成被试 %s，时长 %s。\n\n', subject_number, duration_label);
        end

    catch ME
        fprintf(2, '!!!! 处理被试 %s 时出错: %s\n', subject_number, ME.message);
        fprintf(2, '     错误发生在: %s (行号: %d)\n', ME.stack(1).name, ME.stack(1).line);
        fprintf('跳过被试 %s，继续下一个...\n\n', subject_number);
        continue;
    end
end

fprintf('====== 所有被试和时长处理完毕！======\n');
