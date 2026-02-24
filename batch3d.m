%-----------------------------------------------------------------------
% Job saved on 04-May-2024 21:59:31 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%%

% 自动读取被试号
base_dir = '/Users/alexjiangzhe/Documents/My_project/1level_traning/output/single_trial_model/';
% 获取基础目录下的所有文件夹
dir_contents = dir(base_dir);
% 筛选出文件夹（排除.和..）
folders = dir_contents([dir_contents.isdir] & ~ismember({dir_contents.name}, {'.', '..'}));
% 提取文件夹名称作为被试号
sub_list = {folders.name};

fprintf('找到 %d 个被试: %s\n', length(sub_list), strjoin(sub_list, ', '));

ses_list={'0-3'};% fz0-1.5 fz1.5-3 fz3-4.5 fz4.5-6 ,'1.5-3','3-4.5'

for i_sub=1:length(sub_list)
    sub=sub_list{i_sub};
    fprintf('正在处理被试: %s\n', sub);
    
    for i_ses=1:length(ses_list)
        ses=ses_list{i_ses};
        fprintf('  正在处理会话: %s\n', ses);
        
        file_dir=[base_dir sub '/' ses '/beta_0*.nii'];
        a = dir(file_dir);
        
        % 检查是否找到文件
        if isempty(a)
            fprintf('    警告: 在 %s 中未找到beta文件\n', file_dir);
            continue;
        end
        
        beta_list={};
        % 动态调整文件数量，使用实际找到的文件数
        num_files = min(101, length(a)); % 最多140个，或实际文件数
        fprintf('    找到 %d 个beta文件，将处理前 %d 个\n', length(a), num_files);
        
        for i = 1:num_files
            beta_list{end+1,1} = [a(i).folder '\' a(i).name];
        end
        
        % 清空之前的matlabbatch
        clear matlabbatch;
        
        matlabbatch{1}.spm.util.cat.vols = beta_list;
        matlabbatch{1}.spm.util.cat.name = 'beta.nii';
        matlabbatch{1}.spm.util.cat.dtype = 4;
        matlabbatch{1}.spm.util.cat.RT = NaN;
        
        % 运行批处理
        try
            spm_jobman('run',matlabbatch);
            fprintf('    成功完成被试 %s 会话 %s 的处理\n', sub, ses);
        catch ME
            fprintf('    错误: 处理被试 %s 会话 %s 时出错: %s\n', sub, ses, ME.message);
        end
    end
end

fprintf('所有被试处理完成！\n');