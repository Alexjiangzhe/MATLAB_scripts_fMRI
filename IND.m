% 表征相似性分析 (RSA) 脚本 - 多时间段自动处理版
clear; clc;

% ========== 输入参数设置 ==========
subjects = [1]; % 被试编号
filebase = '/Users/alexjiangzhe/Documents/My_project/1level_traning/output/single_trial_model'; % 基础路径
masksPth = '/Users/alexjiangzhe/Documents/My_project/internship_He_researsh_group'; % ROI基础路径
corr_type = 'Pearson'; % 相似性度量
timeIntervals = {'0-3', '0-1.5', '1.5-3'}; % 添加要处理的时间段

% 通过对话框选择多个ROI文件
[roiFiles, roiPath] = uigetfile(fullfile(masksPth, '*.nii'), ...
    '选择ROI文件(可多选)', 'MultiSelect', 'on');
if isequal(roiFiles, 0)
    error('未选择任何ROI文件');
end
% 确保roiFiles是cell数组
if ischar(roiFiles)
    roisFILES = {roiFiles};
else
    roisFILES = roiFiles;
end

% 实验参数
nTrials = 100; % 总试次数
trialsPerCond = 20; % 每个条件试次数
conditions = {'blank', 'distancing', 'distraction', ...
              'inhibition', ...
              'nature'}; % 条件名称
nConditions = length(conditions);

% 结果根目录
resultRoot = '/Users/alexjiangzhe/Documents/My_project/1level_traning/output/single_trial_model/result';

% ========== 主分析循环 ==========
for tIdx = 1:length(timeIntervals) % 添加时间段循环
    timeInt = timeIntervals{tIdx};
    fprintf('\n===== 处理时间段: %s =====\n', timeInt);
    
    for subIdx = 1:length(subjects)
        sub = subjects(subIdx);
        fprintf('\n处理被试: %d | 时间段: %s\n', sub, timeInt);
        
        % 构建路径
        subPath = fullfile(filebase, num2str(sub)); % G:\deinfo\1st\single_trial_model\6
        betaPath = fullfile(subPath, timeInt); % 使用动态时间段路径
        
        % 创建被试结果目录
        subResultDir = fullfile(resultRoot, num2str(sub));
        if ~exist(subResultDir, 'dir')
            mkdir(subResultDir);
        end
        
        % 创建时间区间目录
        timeIntervalDir = fullfile(subResultDir, timeInt);
        if ~exist(timeIntervalDir, 'dir')
            mkdir(timeIntervalDir);
        end
        
        % 检查beta路径是否存在
        if ~exist(betaPath, 'dir')
            warning('时间段目录不存在: %s，跳过此时间段', betaPath);
            continue;
        end
        
        % 清理之前的con文件
        conFiles = dir(fullfile(betaPath, 'con_*.nii'));
        if ~isempty(conFiles)
            for f = 1:length(conFiles)
                delete(fullfile(conFiles(f).folder, conFiles(f).name));
            end
            fprintf('已删除 %d 个旧con文件\n', length(conFiles));
        end
        
        % 获取有效beta文件 (跳过51-56和147-154)
        betaFiles = {};
        totalScans = 107; % 假设每个时间段有426个扫描
        
        for i = 1:totalScans
            % 跳过不需要的扫描
            if (i >= 101 && i <= 107) % 根据实际情况调整跳过规则
                continue;
            end
            
            % 构建文件名
            fname = sprintf('beta_%04d.nii', i);
            filePath = fullfile(betaPath, fname);
            
            % 检查文件是否存在
            if exist(filePath, 'file')
                betaFiles{end+1} = filePath;
            else
                warning('文件不存在: %s', filePath);
            end
        end
        
        % 验证文件数量
        if length(betaFiles) ~= nTrials
            error('发现 %d 个beta文件, 但需要 %d 个', length(betaFiles), nTrials);
        else
            fprintf('找到 %d 个有效beta文件\n', nTrials);
        end
        
        % 创建con文件 (直接复制beta文件)
        for t = 1:nTrials
            origFile = betaFiles{t};
            [path, name, ext] = fileparts(origFile);
            conFile = fullfile(path, ['con_' name(6:end) ext]); % 从beta_0001创建con_0001
            
            copyfile(origFile, conFile);
        end
        fprintf('已创建 %d 个con文件\n', nTrials);
        
        % ========== 设置RSA分析 ==========
        cfg = decoding_defaults;
        cfg.software = 'SPM12';
        cfg.plot_design = 0;
        cfg.results.overwrite = 1;
        cfg.analysis = 'roi';
        
        % 获取所有con文件
        conFiles = dir(fullfile(betaPath, 'con_*.nii'));
        cfg.files.name = cell(nTrials, 1);
        for t = 1:nTrials
            cfg.files.name{t} = fullfile(conFiles(t).folder, conFiles(t).name);
        end
        
        % 创建标签向量
        cfg.files.label = zeros(nTrials, 1);
        for condIdx = 1:nConditions
            startIdx = (condIdx-1)*trialsPerCond + 1;
            endIdx = condIdx*trialsPerCond;
            cfg.files.label(startIdx:endIdx) = condIdx;
        end
        
        cfg.files.chunk = ones(nTrials, 1); % 所有试次属于同一组块
        
        % 设置相似性分析
        cfg.decoding.software = 'similarity';
        cfg.decoding.method = 'classification';
        cfg.decoding.train.classification.model_parameters = corr_type;
        cfg.results.output = 'other';
        cfg.verbose = 0;
        
        % 创建设计矩阵
        cfg.design = make_design_similarity(cfg);
        cfg.design.unbalanced_data = 'ok';
        
        % ========== ROI批量分析 ==========
        for r = 1:length(roisFILES)
            roiName = roisFILES{r};
            fprintf('\n正在分析ROI (%d/%d): %s\n', r, length(roisFILES), roiName);
            
            % ROI路径
            roiMask = fullfile(roiPath, roiName);
            
            % 检查ROI文件是否存在
            if ~exist(roiMask, 'file')
                warning('ROI文件不存在: %s，跳过此ROI', roiMask);
                continue;
            end
            
            % 创建ROI结果目录
            [~, roiBaseName] = fileparts(roiName); % 去除.nii后缀
            roiResultsDir = fullfile(timeIntervalDir, roiBaseName);
            if ~exist(roiResultsDir, 'dir')
                mkdir(roiResultsDir);
            end
            cfg.results.dir = roiResultsDir;
            cfg.files.mask = roiMask;
            
            % 运行RSA分析
            try
                results = decoding(cfg);
                
                % 处理结果
                similarity_matrix = results.other.output{1};
                rdm = 1 - similarity_matrix; % 表征不相似性矩阵
                
                % 保存完整RDM
                save(fullfile(roiResultsDir, 'full_rdm.mat'), 'rdm');
                fprintf('完整RDM已保存至: %s\n', roiResultsDir);
                
                % ===== 为每个条件创建并保存35x35的RDM =====
                for condIdx = 1:nConditions
                    % 获取当前条件的索引范围
                    startIdx = (condIdx-1)*trialsPerCond + 1;
                    endIdx = condIdx*trialsPerCond;
                    condIndices = startIdx:endIdx;
                    
                    % 提取当前条件的35x35子矩阵
                    cond_rdm = rdm(condIndices, condIndices);
                    
                    % 保存条件RDM
                    condName = conditions{condIdx};
                    % 替换条件名称中的特殊字符为下划线
                    safeCondName = regexprep(condName, '[^\w]', '_');
                    condFileName = fullfile(roiResultsDir, ['rdm_', safeCondName, '.mat']);
                    
                    save(condFileName, 'cond_rdm');
                    fprintf('条件 "%s" 的35x35 RDM已保存: %s\n', condName, condFileName);
                    
                    % ===== 可视化当前条件的RDM =====
                    figure('Position', [100, 100, 600, 550], 'Visible', 'off');
                    imagesc(cond_rdm);
                    colormap('jet');
                    colorbar;
                    
                    % 设置标题
                    title(sprintf('Sub-%d | %s | %s | 条件: %s', sub, timeInt, roiBaseName, condName), 'FontSize', 14);
                    
                    % 添加轴标签（试次编号）
                    ax = gca;
                    ax.XTick = 1:trialsPerCond;
                    ax.YTick = 1:trialsPerCond;
                    ax.FontSize = 8;
                    axis square;
                    
                    % 添加网格线
                    hold on;
                    for k = 1:trialsPerCond
                        plot([k+0.5, k+0.5], [0.5, trialsPerCond+0.5], 'k-', 'LineWidth', 0.5);
                        plot([0.5, trialsPerCond+0.5], [k+0.5, k+0.5], 'k-', 'LineWidth', 0.5);
                    end
                    
                    % 保存图像
                    condImageName = fullfile(roiResultsDir, ['rdm_', safeCondName, '.jpg']);
                    saveas(gcf, condImageName);
                    close(gcf);
                end
                
                % ===== 可视化完整RDM =====
                figure('Position', [100, 100, 800, 700], 'Visible', 'off');
                imagesc(rdm);
                colormap('jet');
                colorbar;
                title(sprintf('Sub-%d | %s | %s | 完整RDM', sub, timeInt, roiBaseName), 'FontSize', 14);
                
                % 添加条件分隔线
                hold on;
                for k = 1:(nConditions-1)
                    line_pos = k*trialsPerCond + 0.5;
                    plot([line_pos, line_pos], [0.5, nTrials+0.5], 'k-', 'LineWidth', 1.5);
                    plot([0.5, nTrials+0.5], [line_pos, line_pos], 'k-', 'LineWidth', 1.5);
                end
                
                % 添加标签
                ax = gca;
                ax.XTick = trialsPerCond/2 + trialsPerCond*(0:(nConditions-1));
                ax.YTick = trialsPerCond/2 + trialsPerCond*(0:(nConditions-1));
                ax.XTickLabel = conditions;
                ax.YTickLabel = conditions;
                ax.FontSize = 10;
                axis square;
                
                % 保存图像
                saveas(gcf, fullfile(roiResultsDir, 'full_rdm_matrix.jpg'));
                saveas(gcf, fullfile(roiResultsDir, 'full_rdm_matrix.fig'));
                close(gcf);
                
            catch ME
                warning('ROI分析失败: %s\n错误信息: %s', roiName, ME.message);
                continue;
            end
        end
    end
end
fprintf('\n所有分析完成!\n');