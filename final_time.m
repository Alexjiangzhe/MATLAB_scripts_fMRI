% 修正后的fMRI数据提取代码（分文件夹存储，独立编号和CSV）
clear; clc;
spm('defaults', 'FMRI');
spm_jobman('initcfg');

%% 路径设置
data_base      = 'D:\unsmooth\combine';
onset_file     = 'D:\unsmooth\combine\onsets_audi_aligned_by_condition.mat';
roi_template   = 'D:\MATLAB21\toolbox\spm12\toolbox\AAL3\ROI_MNI_V7.nii';
roi_label_xml  = 'D:\MATLAB21\toolbox\spm12\toolbox\AAL3\ROI_MNI_V7.xml';
roi_14net_path = 'D:\unsmooth\14net';  % 14net ROI路径
output_path    = 'D:\unsmooth\combine\result';

% 创建子文件夹
wholebrain_output_path = fullfile(output_path, 'WholeBrain');
aal_output_path = fullfile(output_path, 'AAL_ROI');
net14_output_path = fullfile(output_path, '14net_ROI');
excel_output_path = fullfile(output_path, 'excel_subjects');  % Excel输出路径

% 创建所有必要的文件夹
folders_to_create = {output_path, wholebrain_output_path, aal_output_path, net14_output_path, excel_output_path};
for i = 1:length(folders_to_create)
    if ~exist(folders_to_create{i}, 'dir')
        mkdir(folders_to_create{i});
        fprintf('创建文件夹: %s\n', folders_to_create{i});
    end
end

%% 严格定义要处理的被试ID列表
subject_ids = {'20','1','2','5','6','7','8','9','13','14','15','16','17','18','19','21','22','23','24','25','26','27','28','29','30','31','32'};

%% 检查文件是否存在
fprintf('=== 文件检查 ===\n');
files_to_check = {roi_template, roi_label_xml, onset_file, roi_14net_path};
file_names = {'AAL ROI模板', 'AAL ROI标签', 'Onset文件', '14net ROI路径'};
for i = 1:length(files_to_check)
    if exist(files_to_check{i}, 'file') || exist(files_to_check{i}, 'dir')
        fprintf('✓ %s存在: %s\n', file_names{i}, files_to_check{i});
    else
        error('✗ %s不存在: %s', file_names{i}, files_to_check{i});
    end
end

%% 加载14net ROI文件
fprintf('\n=== 加载14net ROI文件 ===\n');
roi_14net_files = dir(fullfile(roi_14net_path, '*.nii'));
if isempty(roi_14net_files)
    error('14net路径下没有找到nii文件: %s', roi_14net_path);
end

roi_14net_info = struct();
for i = 1:length(roi_14net_files)
    [~, roi_name, ~] = fileparts(roi_14net_files(i).name);
    roi_name = matlab.lang.makeValidName(roi_name);  % 确保变量名合法
    roi_14net_info(i).name = roi_name;
    roi_14net_info(i).file = fullfile(roi_14net_files(i).folder, roi_14net_files(i).name);
    fprintf('发现14net ROI %d: %s\n', i, roi_name);
end
fprintf('总共找到 %d 个14net ROI文件\n', length(roi_14net_info));

%% 加载触发时间 - 修复版本
fprintf('\n=== 加载 onsets 文件 ===\n');
try
    % 首先检查文件中的变量
    file_vars = whos('-file', onset_file);
    var_names = {file_vars.name};
    fprintf('onset文件中的变量: %s\n', strjoin(var_names, ', '));
    
    % 加载文件
    loaded_data = load(onset_file);
    
    % 智能查找onsets变量
    onsets_var_name = '';
    if ismember('onsets_audi_aligned_by_condition', var_names)
        onsets_var_name = 'onsets_audi_aligned_by_condition';
        onsets_audi_aligned_by_condition = loaded_data.onsets_audi_aligned_by_condition;
    elseif ismember('onsets_auditory', var_names)
        onsets_var_name = 'onsets_auditory';
        onsets_audi_aligned_by_condition = loaded_data.onsets_auditory;  % 重命名为统一变量
    elseif ismember('onsets_visual', var_names)
        onsets_var_name = 'onsets_visual';
        onsets_audi_aligned_by_condition = loaded_data.onsets_visual;  % 重命名为统一变量
    else
        % 尝试找到包含'onset'的变量
        for i = 1:length(var_names)
            if contains(lower(var_names{i}), 'onset')
                onsets_var_name = var_names{i};
                onsets_audi_aligned_by_condition = loaded_data.(var_names{i});
                break;
            end
        end
    end
    
    if isempty(onsets_var_name)
        error('未找到onsets变量，请检查文件内容');
    else
        fprintf('使用变量: %s\n', onsets_var_name);
    end
    
catch ME
    error('onsets文件加载失败: %s', ME.message);
end

% 检查onsets文件中的被试是否与subject_ids匹配
if ~isstruct(onsets_audi_aligned_by_condition)
    error('onsets变量不是结构体格式');
end

all_subs_in_onsets = fieldnames(onsets_audi_aligned_by_condition);
fprintf('onsets文件中包含的被试总数: %d\n', length(all_subs_in_onsets));

% 找出需要处理的被试（存在于subject_ids中的）
valid_subs = {};
for i = 1:length(subject_ids)
    sub_str = ['sub' subject_ids{i}];
    if ismember(sub_str, all_subs_in_onsets)
        valid_subs{end+1} = subject_ids{i};
    else
        warning('被试 %s 在onsets文件中不存在', subject_ids{i});
    end
end

fprintf('实际需要处理的被试: %s\n', strjoin(valid_subs, ', '));
if isempty(valid_subs)
    error('没有找到需要处理的被试！');
end

%% 加载 AAL ROI 模板和标签
fprintf('\n=== 加载AAL ROI数据 ===\n');
try
    % 加载ROI模板
    roiV = spm_vol(roi_template);
    roi_data = spm_read_vols(roiV);
    fprintf('AAL ROI模板维度: [%s]\n', num2str(size(roi_data)));
    
    % 解析XML标签文件
    labelDoc = xmlread(roi_label_xml);
    labelItems = labelDoc.getElementsByTagName('label');
    
    aal_roi_name_map = containers.Map('KeyType','double','ValueType','char');
    for i = 0:labelItems.getLength-1
        item = labelItems.item(i);
        id_node = item.getElementsByTagName('index').item(0);
        id_val = str2double(char(id_node.getFirstChild.getNodeValue));
        name_node = item.getElementsByTagName('name').item(0);
        roi_name = char(name_node.getFirstChild.getNodeValue);
        aal_roi_name_map(id_val) = strtrim(roi_name);
    end
    aal_roi_id_list = cell2mat(aal_roi_name_map.keys);
    fprintf('成功加载 %d 个AAL ROI区域\n', length(aal_roi_id_list));
catch ME
    error('AAL ROI数据加载失败: %s', ME.message);
end

%% 创建分类的ROI对照表
fprintf('\n=== 创建分类对照表 ===\n');

% 收集所有条件名称（从第一个被试获取）
first_sub_str = ['sub' valid_subs{1}];
all_conditions = fieldnames(onsets_audi_aligned_by_condition.(first_sub_str));

% 创建condition对照表（所有类型共用）
condition_mapping = {};
for i = 1:length(all_conditions)
    condition_mapping{i, 1} = i;
    condition_mapping{i, 2} = all_conditions{i};
end

% 创建全脑ROI对照表（编号从1开始）
wholebrain_roi_mapping = {};
wholebrain_roi_mapping{1, 1} = 1;
wholebrain_roi_mapping{1, 2} = 'WholeBrain_Average';

% 创建AAL ROI对照表（编号从1开始）
aal_roi_mapping = {};
counter = 1;
for roi_id = aal_roi_id_list
    aal_roi_mapping{counter, 1} = counter;
    aal_roi_mapping{counter, 2} = aal_roi_name_map(roi_id);
    counter = counter + 1;
end

% 创建14net ROI对照表（编号从1开始）
net14_roi_mapping = {};
for i = 1:length(roi_14net_info)
    net14_roi_mapping{i, 1} = i;
    net14_roi_mapping{i, 2} = roi_14net_info(i).name;
end

% 保存各类型的对照表
try
    % 条件对照表保存到总文件夹
    condition_table = cell2table(condition_mapping, 'VariableNames', {'条件编号', '条件名称'});
    condition_lookup_path = fullfile(excel_output_path, 'Condition_对照表.xlsx');
    writetable(condition_table, condition_lookup_path, 'Sheet', 'Condition对照表', 'WriteRowNames', false);
    
    % 全脑ROI对照表
    wb_roi_table = cell2table(wholebrain_roi_mapping, 'VariableNames', {'ROI编号', 'ROI名称'});
    wb_roi_lookup_path = fullfile(wholebrain_output_path, 'WholeBrain_ROI_对照表.xlsx');
    writetable(wb_roi_table, wb_roi_lookup_path, 'Sheet', 'ROI对照表', 'WriteRowNames', false);
    
    % AAL ROI对照表
    aal_roi_table = cell2table(aal_roi_mapping, 'VariableNames', {'ROI编号', 'ROI名称'});
    aal_roi_lookup_path = fullfile(aal_output_path, 'AAL_ROI_对照表.xlsx');
    writetable(aal_roi_table, aal_roi_lookup_path, 'Sheet', 'ROI对照表', 'WriteRowNames', false);
    
    % 14net ROI对照表
    net14_roi_table = cell2table(net14_roi_mapping, 'VariableNames', {'ROI编号', 'ROI名称'});
    net14_roi_lookup_path = fullfile(net14_output_path, '14net_ROI_对照表.xlsx');
    writetable(net14_roi_table, net14_roi_lookup_path, 'Sheet', 'ROI对照表', 'WriteRowNames', false);
    
    fprintf('✓ 分类对照表创建成功\n');
    fprintf('  全脑ROI数量: %d\n', size(wholebrain_roi_mapping, 1));
    fprintf('  AAL ROI数量: %d\n', size(aal_roi_mapping, 1));
    fprintf('  14net ROI数量: %d\n', size(net14_roi_mapping, 1));
    fprintf('  实验条件数量: %d\n', length(all_conditions));
    
catch ME
    warning('对照表创建失败: %s', ME.message);
end

%% 参数设置
TR = 1.5;
n_valid_scans = 4;

%% 初始化CSV数据存储（分类存储）
csv_headers = {'Subject', 'Condition', 'Trial', 'ROI_Name', 'TimePoint', 'Signal'};
wholebrain_csv_data = {};  % 全脑CSV数据
aal_csv_data = {};         % AAL CSV数据
net14_csv_data = {};       % 14net CSV数据

%% 主处理循环
fprintf('\n===== 开始处理被试数据 =====\n');
fprintf('总被试数: %d\n', length(valid_subs));

for sub_idx = 1:length(valid_subs)
    subject_id = valid_subs{sub_idx};
    sub_str = ['sub' subject_id];
    
    fprintf('\n=== 处理进度: %d/%d ===\n', sub_idx, length(valid_subs));
    fprintf('=== 处理被试 %s ===\n', subject_id);
    tic;
    
    %% 初始化当前被试的Excel数据存储
    wholebrain_excel_data = [];
    aal_excel_data = [];
    net14_excel_data = [];
    
    %% 检查被试数据
    sub_folder = fullfile(data_base, subject_id);
    nii_file = fullfile(sub_folder, 'merged_4d.nii');
    if ~exist(nii_file, 'file')
        warning('NIfTI文件缺失: %s，跳过', nii_file);
        continue;
    end
    
    %% 加载fMRI数据
    try
        V = spm_vol(nii_file);
        dims = V(1).dim;
        num_volumes = numel(V);
        fprintf('fMRI数据维度: [%s]，时间点: %d\n', num2str(dims), num_volumes);
    catch ME
        warning('fMRI数据加载失败: %s，跳过', ME.message);
        continue;
    end
    
    %% 检查维度匹配
    aal_roi_compatible = isequal(size(roi_data), dims);
    if ~aal_roi_compatible
        warning('AAL ROI/fMRI维度不匹配，跳过AAL ROI分析');
    end
    
    %% 加载14net ROI数据
    roi_14net_data = {};
    roi_14net_compatible = false(length(roi_14net_info), 1);
    
    fprintf('加载14net ROI数据...\n');
    for i = 1:length(roi_14net_info)
        try
            roi_vol = spm_vol(roi_14net_info(i).file);
            roi_mask = spm_read_vols(roi_vol);
            if isequal(size(roi_mask), dims)
                roi_14net_data{i} = roi_mask;
                roi_14net_compatible(i) = true;
                fprintf('  ✓ %s 加载成功\n', roi_14net_info(i).name);
            else
                warning('  ✗ %s 维度不匹配，跳过', roi_14net_info(i).name);
            end
        catch ME
            warning('  ✗ %s 加载失败: %s', roi_14net_info(i).name, ME.message);
        end
    end
    
    %% 获取当前被试的条件信息
    if ~isfield(onsets_audi_aligned_by_condition, sub_str)
        warning('被试 %s 在onsets中不存在，跳过', subject_id);
        continue;
    end
    
    condition_labels = fieldnames(onsets_audi_aligned_by_condition.(sub_str));
    fprintf('实验条件 (%d): %s\n', length(condition_labels), strjoin(condition_labels, ', '));
    
    %% 初始化存储结构
    voxel_result = struct();
    aal_roi_result = struct();
    roi_14net_result = struct();
    
    for c = 1:length(condition_labels)
        cond = matlab.lang.makeValidName(condition_labels{c});
        voxel_result.(cond) = [];
        
        if aal_roi_compatible
            for i = 1:size(aal_roi_mapping, 1)
                roi_name = matlab.lang.makeValidName(aal_roi_mapping{i, 2});
                aal_roi_result.(roi_name).(cond) = [];
            end
        end
        
        for i = 1:length(roi_14net_info)
            if roi_14net_compatible(i)
                roi_name = roi_14net_info(i).name;
                roi_14net_result.(roi_name).(cond) = [];
            end
        end
    end
    
    %% 处理每个条件
    for cond_idx = 1:length(condition_labels)
        cond_name = condition_labels{cond_idx};
        cond_field = matlab.lang.makeValidName(cond_name);
        onsets_original = onsets_audi_aligned_by_condition.(sub_str).(cond_name);
        
        % 确保数据是数值型
        if iscell(onsets_original)
            onsets_original = cell2mat(onsets_original);
        end
        
        % 去除NaN值和无效值
        valid_onset_mask = ~isnan(onsets_original) & onsets_original > 0;
        onsets_clean = onsets_original(valid_onset_mask);
        
        fprintf('\n条件 %d/%d: %s (原始:%d trials, 有效:%d trials)\n', ...
               cond_idx, length(condition_labels), cond_name, length(onsets_original), length(onsets_clean));
        
        %% 检查哪些trial的时间窗口是有效的
        valid_trials = 0;
        vol_data = zeros([dims, n_valid_scans, length(onsets_clean)], 'single');
        valid_trial_mapping = [];  % 存储有效trial在原始数组中的索引
        valid_onsets_final = [];   % 存储最终有效的onset值
        
        for t = 1:length(onsets_clean)
            scan_idx = floor(onsets_clean(t)/TR) + (1:n_valid_scans);
            if any(scan_idx > num_volumes) || any(scan_idx < 1)
                fprintf('  Trial %d 跳过 (onset=%.2f, 扫描索引超出范围)\n', t, onsets_clean(t));
                continue;
            end
            
            valid_trials = valid_trials + 1;
            valid_trial_mapping(end+1) = t;  % 记录在原始clean数组中的位置
            valid_onsets_final(end+1) = onsets_clean(t);  % 记录对应的onset值
            
            for tt = 1:n_valid_scans
                vol_data(:,:,:,tt,valid_trials) = single(spm_read_vols(V(scan_idx(tt))));
            end
        end
        
        fprintf('  最终有效trial数: %d\n', valid_trials);
        if valid_trials > 0
            fprintf('  有效onset值: ');
            for i = 1:min(5, length(valid_onsets_final))
                fprintf('%.2f ', valid_onsets_final(i));
            end
            fprintf('\n');
        end
        
        if valid_trials > 0
            voxel_result.(cond_field) = vol_data(:,:,:,:,1:valid_trials);
            fprintf('  保存体素数据: [%s]\n', num2str(size(voxel_result.(cond_field))));
            
            % 提取全脑平均信号并添加到Excel数据和CSV数据
            for t = 1:valid_trials
                trial_idx_in_original = valid_trial_mapping(t);  % 在原始数组中的位置
                onset_value = valid_onsets_final(t);  % 对应的真实onset值
                
                fprintf('    处理Trial %d (原始位置:%d, onset:%.2f)\n', t, trial_idx_in_original, onset_value);
                
                % 获取4个时间点的全脑信号
                brain_signals = zeros(1, 4);
                for tt = 1:n_valid_scans
                    brain_signal = mean(vol_data(:,:,:,tt,t), 'all', 'omitnan');
                    brain_signals(tt) = brain_signal;
                    
                    % 添加到全脑CSV数据
                    wholebrain_csv_data{end+1,1} = subject_id;
                    wholebrain_csv_data{end,2} = cond_name;
                    wholebrain_csv_data{end,3} = trial_idx_in_original;
                    wholebrain_csv_data{end,4} = 'WholeBrain_Average';
                    wholebrain_csv_data{end,5} = tt;
                    wholebrain_csv_data{end,6} = brain_signal;
                end
                
                % 添加到全脑Excel数据（使用新的编号系统）
                network_num = map_to_number('WholeBrain_Average', wholebrain_roi_mapping);
                condition_num = map_to_number(cond_name, condition_mapping);
                subject_num = str2double(subject_id);
                
                row_data = [network_num, condition_num, subject_num, trial_idx_in_original, ...
                           brain_signals(1), brain_signals(2), brain_signals(3), brain_signals(4), onset_value];
                wholebrain_excel_data = [wholebrain_excel_data; row_data];
                
                fprintf('      添加到全脑Excel: onset=%.2f\n', onset_value);
            end
        end
        
        %% AAL ROI数据提取
        if aal_roi_compatible && valid_trials > 0
            fprintf('  提取AAL ROI信号...');
            counter = 1;
            for roi_id = aal_roi_id_list
                mask = (roi_data == roi_id);
                if nnz(mask) == 0
                    counter = counter + 1;
                    continue;
                end
                
                roi_name_original = aal_roi_name_map(roi_id);
                roi_name = matlab.lang.makeValidName(roi_name_original);
                roi_ts = zeros(valid_trials, n_valid_scans);
                
                for t = 1:valid_trials
                    trial_idx_in_original = valid_trial_mapping(t);
                    onset_value = valid_onsets_final(t);
                    roi_signals = zeros(1, 4);
                    
                    for tt = 1:n_valid_scans
                        tmp = vol_data(:,:,:,tt,t);
                        roi_signal = mean(tmp(mask), 'omitnan');
                        roi_ts(t,tt) = roi_signal;
                        roi_signals(tt) = roi_signal;
                        
                        % 添加到AAL CSV数据
                        aal_csv_data{end+1,1} = subject_id;
                        aal_csv_data{end,2} = cond_name;
                        aal_csv_data{end,3} = trial_idx_in_original;
                        aal_csv_data{end,4} = roi_name_original;
                        aal_csv_data{end,5} = tt;
                        aal_csv_data{end,6} = roi_signal;
                    end
                    
                    % 添加到AAL Excel数据（使用新的编号系统）
                    network_num = counter;  % 直接使用counter作为编号
                    condition_num = map_to_number(cond_name, condition_mapping);
                    subject_num = str2double(subject_id);
                    
                    row_data = [network_num, condition_num, subject_num, trial_idx_in_original, ...
                               roi_signals(1), roi_signals(2), roi_signals(3), roi_signals(4), onset_value];
                    aal_excel_data = [aal_excel_data; row_data];
                end
                
                aal_roi_result.(roi_name).(cond_field) = roi_ts;
                counter = counter + 1;
            end
            fprintf('完成\n');
        end
        
        %% 14net ROI数据提取
        if valid_trials > 0
            fprintf('  提取14net ROI信号...');
            for i = 1:length(roi_14net_info)
                if ~roi_14net_compatible(i), continue; end
                
                roi_name = roi_14net_info(i).name;
                mask = roi_14net_data{i} > 0;  % 假设非零值为ROI区域
                if nnz(mask) == 0, continue; end
                
                roi_ts = zeros(valid_trials, n_valid_scans);
                
                for t = 1:valid_trials
                    trial_idx_in_original = valid_trial_mapping(t);
                    onset_value = valid_onsets_final(t);
                    roi_signals = zeros(1, 4);
                    
                    for tt = 1:n_valid_scans
                        tmp = vol_data(:,:,:,tt,t);
                        roi_signal = mean(tmp(mask), 'omitnan');
                        roi_ts(t,tt) = roi_signal;
                        roi_signals(tt) = roi_signal;
                        
                        % 添加到14net CSV数据
                        net14_csv_data{end+1,1} = subject_id;
                        net14_csv_data{end,2} = cond_name;
                        net14_csv_data{end,3} = trial_idx_in_original;
                        net14_csv_data{end,4} = roi_name;
                        net14_csv_data{end,5} = tt;
                        net14_csv_data{end,6} = roi_signal;
                    end
                    
                    % 添加到14net Excel数据（使用新的编号系统）
                    network_num = map_to_number(roi_name, net14_roi_mapping);
                    condition_num = map_to_number(cond_name, condition_mapping);
                    subject_num = str2double(subject_id);
                    
                    row_data = [network_num, condition_num, subject_num, trial_idx_in_original, ...
                               roi_signals(1), roi_signals(2), roi_signals(3), roi_signals(4), onset_value];
                    net14_excel_data = [net14_excel_data; row_data];
                end
                
                roi_14net_result.(roi_name).(cond_field) = roi_ts;
            end
            fprintf('完成\n');
        end
    end
    
    %% 保存当前被试的Excel文件到不同文件夹
    fprintf('\n保存Excel文件到分类文件夹...\n');
    
    % 保存全脑Excel文件
    if ~isempty(wholebrain_excel_data)
        excel_filename = sprintf('%s.xlsx', subject_id);
        excel_filepath = fullfile(wholebrain_output_path, excel_filename);
        
        try
            xlswrite(excel_filepath, wholebrain_excel_data);
            fprintf('✓ 全脑Excel文件: %s (数据: %d行)\n', excel_filename, size(wholebrain_excel_data, 1));
        catch ME
            warning('全脑Excel文件创建失败: %s', ME.message);
        end
    end
    
    % 保存AAL Excel文件
    if ~isempty(aal_excel_data)
        excel_filename = sprintf('%s.xlsx', subject_id);
        excel_filepath = fullfile(aal_output_path, excel_filename);
        
        try
            xlswrite(excel_filepath, aal_excel_data);
            fprintf('✓ AAL Excel文件: %s (数据: %d行)\n', excel_filename, size(aal_excel_data, 1));
        catch ME
            warning('AAL Excel文件创建失败: %s', ME.message);
        end
    end
    
    % 保存14net Excel文件
    if ~isempty(net14_excel_data)
        excel_filename = sprintf('%s.xlsx', subject_id);
        excel_filepath = fullfile(net14_output_path, excel_filename);
        
        try
            xlswrite(excel_filepath, net14_excel_data);
            fprintf('✓ 14net Excel文件: %s (数据: %d行)\n', excel_filename, size(net14_excel_data, 1));
        catch ME
            warning('14net Excel文件创建失败: %s', ME.message);
        end
    end
    
    %% 保存.mat结果文件到相应文件夹
    fprintf('\n保存.mat结果文件到分类文件夹...\n');
    try
        % 保存体素数据到全脑文件夹
        voxel_file = fullfile(wholebrain_output_path, sprintf('voxel_sub%s.mat', subject_id));
        save(voxel_file, 'voxel_result', '-v7.3');
        fprintf('✓ 体素数据保存至: %s\n', voxel_file);
        
        % 保存AAL ROI数据到AAL文件夹
        if aal_roi_compatible
            aal_roi_file = fullfile(aal_output_path, sprintf('aal_roi_sub%s.mat', subject_id));
            save(aal_roi_file, 'aal_roi_result', '-v7.3');
            fprintf('✓ AAL ROI数据保存至: %s\n', aal_roi_file);
        end
        
        % 保存14net ROI数据到14net文件夹
        if any(roi_14net_compatible)
            roi_14net_file = fullfile(net14_output_path, sprintf('14net_roi_sub%s.mat', subject_id));
            save(roi_14net_file, 'roi_14net_result', '-v7.3');
            fprintf('✓ 14net ROI数据保存至: %s\n', roi_14net_file);
        end
    catch ME
        warning('保存.mat文件失败: %s', ME.message);
    end
    
    fprintf('被试处理完成! 耗时: %.2f 秒\n', toc);
end

%% 保存分类CSV文件
fprintf('\n=== 保存分类CSV文件 ===\n');

% 保存全脑CSV文件
try
    if ~isempty(wholebrain_csv_data)
        csv_file = fullfile(wholebrain_output_path, 'WholeBrain_timeseries.csv');
        write_csv_file(csv_file, csv_headers, wholebrain_csv_data);
        fprintf('✓ 全脑CSV文件保存成功: %s (数据行数: %d)\n', csv_file, size(wholebrain_csv_data, 1));
    end
catch ME
    warning('全脑CSV文件保存失败: %s', ME.message);
end

% 保存AAL CSV文件
try
    if ~isempty(aal_csv_data)
        csv_file = fullfile(aal_output_path, 'AAL_ROI_timeseries.csv');
        write_csv_file(csv_file, csv_headers, aal_csv_data);
        fprintf('✓ AAL CSV文件保存成功: %s (数据行数: %d)\n', csv_file, size(aal_csv_data, 1));
    end
catch ME
    warning('AAL CSV文件保存失败: %s', ME.message);
end

% 保存14net CSV文件
try
    if ~isempty(net14_csv_data)
        csv_file = fullfile(net14_output_path, '14net_ROI_timeseries.csv');
        write_csv_file(csv_file, csv_headers, net14_csv_data);
        fprintf('✓ 14net CSV文件保存成功: %s (数据行数: %d)\n', csv_file, size(net14_csv_data, 1));
    end
catch ME
    warning('14net CSV文件保存失败: %s', ME.message);
end

%% 合并所有CSV数据并保存总文件
fprintf('\n=== 保存合并CSV文件 ===\n');
try
    % 合并所有数据，添加ROI类型列
    all_csv_data = {};
    combined_headers = {'Subject', 'Condition', 'Trial', 'ROI_Type', 'ROI_Name', 'TimePoint', 'Signal'};
    
    % 添加全脑数据
    for i = 1:size(wholebrain_csv_data, 1)
        all_csv_data{end+1,1} = wholebrain_csv_data{i,1};  % Subject
        all_csv_data{end,2} = wholebrain_csv_data{i,2};    % Condition
        all_csv_data{end,3} = wholebrain_csv_data{i,3};    % Trial
        all_csv_data{end,4} = 'WholeBrain';                % ROI_Type
        all_csv_data{end,5} = wholebrain_csv_data{i,4};    % ROI_Name
        all_csv_data{end,6} = wholebrain_csv_data{i,5};    % TimePoint
        all_csv_data{end,7} = wholebrain_csv_data{i,6};    % Signal
    end
    
    % 添加AAL数据
    for i = 1:size(aal_csv_data, 1)
        all_csv_data{end+1,1} = aal_csv_data{i,1};         % Subject
        all_csv_data{end,2} = aal_csv_data{i,2};           % Condition
        all_csv_data{end,3} = aal_csv_data{i,3};           % Trial
        all_csv_data{end,4} = 'AAL_ROI';                   % ROI_Type
        all_csv_data{end,5} = aal_csv_data{i,4};           % ROI_Name
        all_csv_data{end,6} = aal_csv_data{i,5};           % TimePoint
        all_csv_data{end,7} = aal_csv_data{i,6};           % Signal
    end
    
    % 添加14net数据
    for i = 1:size(net14_csv_data, 1)
        all_csv_data{end+1,1} = net14_csv_data{i,1};       % Subject
        all_csv_data{end,2} = net14_csv_data{i,2};         % Condition
        all_csv_data{end,3} = net14_csv_data{i,3};         % Trial
        all_csv_data{end,4} = '14net_ROI';                 % ROI_Type
        all_csv_data{end,5} = net14_csv_data{i,4};         % ROI_Name
        all_csv_data{end,6} = net14_csv_data{i,5};         % TimePoint
        all_csv_data{end,7} = net14_csv_data{i,6};         % Signal
    end
    
    % 保存合并的CSV文件
    combined_csv_file = fullfile(output_path, 'All_Subjects_Combined_timeseries.csv');
    write_csv_file(combined_csv_file, combined_headers, all_csv_data);
    fprintf('✓ 合并CSV文件保存成功: %s (总数据行数: %d)\n', combined_csv_file, size(all_csv_data, 1));
    
catch ME
    warning('合并CSV文件保存失败: %s', ME.message);
end

%% 生成数据摘要
fprintf('\n=== 数据处理摘要 ===\n');
fprintf('成功处理的被试数: %d\n', length(valid_subs));
fprintf('数据统计:\n');
fprintf('  全脑数据行数: %d\n', size(wholebrain_csv_data, 1));
fprintf('  AAL ROI数据行数: %d\n', size(aal_csv_data, 1));
fprintf('  14net ROI数据行数: %d\n', size(net14_csv_data, 1));
fprintf('  总合并数据行数: %d\n', size(all_csv_data, 1));

%% 文件夹结构和文件列表
fprintf('\n=== 文件夹结构 ===\n');
fprintf('主结果文件夹: %s\n', output_path);
fprintf('├── WholeBrain/     (全脑数据，ROI编号1)\n');
fprintf('├── AAL_ROI/        (AAL ROI数据，ROI编号1-%d)\n', size(aal_roi_mapping, 1));
fprintf('├── 14net_ROI/      (14net ROI数据，ROI编号1-%d)\n', size(net14_roi_mapping, 1));
fprintf('├── excel_subjects/ (条件对照表)\n');
fprintf('└── All_Subjects_Combined_timeseries.csv (合并数据)\n');

% 检查各文件夹中的文件
fprintf('\n=== 生成的文件清单 ===\n');

% 全脑文件夹
fprintf('\n全脑文件夹 (%s):\n', wholebrain_output_path);
wb_files = dir(fullfile(wholebrain_output_path, '*'));
wb_files = wb_files(~[wb_files.isdir]);
for i = 1:length(wb_files)
    fprintf('  %s\n', wb_files(i).name);
end

% AAL文件夹
fprintf('\nAAL ROI文件夹 (%s):\n', aal_output_path);
aal_files = dir(fullfile(aal_output_path, '*'));
aal_files = aal_files(~[aal_files.isdir]);
for i = 1:length(aal_files)
    fprintf('  %s\n', aal_files(i).name);
end

% 14net文件夹
fprintf('\n14net ROI文件夹 (%s):\n', net14_output_path);
net14_files = dir(fullfile(net14_output_path, '*'));
net14_files = net14_files(~[net14_files.isdir]);
for i = 1:length(net14_files)
    fprintf('  %s\n', net14_files(i).name);
end

% Excel对照表文件夹
fprintf('\nExcel对照表文件夹 (%s):\n', excel_output_path);
excel_files = dir(fullfile(excel_output_path, '*.xlsx'));
for i = 1:length(excel_files)
    fprintf('  %s\n', excel_files(i).name);
end

fprintf('\n=== 文件说明 ===\n');
fprintf('1. 文件夹分类和编号系统:\n');
fprintf('   - WholeBrain/: 全脑平均信号，ROI编号=1\n');
fprintf('   - AAL_ROI/: AAL图谱ROI，ROI编号=1-%d\n', size(aal_roi_mapping, 1));
fprintf('   - 14net_ROI/: 14网络ROI，ROI编号=1-%d\n', size(net14_roi_mapping, 1));
fprintf('   - excel_subjects/: 条件对照表\n');
fprintf('\n2. 每个子文件夹包含:\n');
fprintf('   - [被试号].xlsx: Excel数据文件\n');
fprintf('   - [类型]_timeseries.csv: CSV数据文件\n');
fprintf('   - [类型]_ROI_对照表.xlsx: ROI编号对照表\n');
fprintf('   - 相应的.mat文件\n');
fprintf('\n3. Excel文件列顺序: network, condition, subject, trial, tr1, tr2, tr3, tr4, onsets\n');
fprintf('   - network列: 在各子文件夹中从1开始重新编号\n');
fprintf('   - onsets列: 包含真实的onset时间值\n');
fprintf('   - 内容: 全部为数字，无表头\n');
fprintf('\n4. CSV文件列顺序: Subject, Condition, Trial, ROI_Name, TimePoint, Signal\n');
fprintf('   - 各子文件夹独立CSV + 总合并CSV\n');
fprintf('   - 合并CSV增加ROI_Type列区分数据类型\n');

fprintf('\n===== 所有被试处理完成 =====\n');
fprintf('结果文件位置: %s\n', output_path);
fprintf('请检查各子文件夹中的文件和对照表是否正确生成。\n');

% 嵌套函数定义
function num_code = map_to_number(value, mapping_cell)
    if iscell(value)
        value = value{1};
    end
    value = char(value);
    
    for i = 1:size(mapping_cell, 1)
        if strcmp(char(mapping_cell{i, 2}), value)
            num_code = mapping_cell{i, 1};
            return;
        end
    end
    num_code = NaN;  % 如果找不到匹配，返回NaN
end

function write_csv_file(filename, headers, data)
    % 写入CSV文件的辅助函数
    fid = fopen(filename, 'w');
    if fid == -1
        error('无法创建CSV文件: %s', filename);
    end
    
    % 写入表头
    fprintf(fid, '%s\n', strjoin(headers, ','));
    
    % 写入数据
    for i = 1:size(data, 1)
        fprintf(fid, '%s,%s,%d,%s,%d,%.6f\n', ...
            data{i,1}, data{i,2}, data{i,3}, data{i,4}, data{i,5}, data{i,6});
    end
    
    fclose(fid);
end