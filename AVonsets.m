% 创建结构来保存onsets和durations
durations_auditory = struct();

% 固定持续时间为4.5秒
fixed_duration = 4.5;

% 听觉组的被试列表
auditory_subjects = {'sub1'};
%                  , 'sub2', 'sub3', 'sub4', 'sub5', 'sub6', 'sub7', 'sub8', ...
%                      'sub9', 'sub10', 'sub11', 'sub12', 'sub13', 'sub14', 'sub15', 'sub16', ...
%                      'sub17', 'sub18', 'sub19', 'sub20', 'sub21', 'sub22', 'sub23', 'sub24', ...
%                      'sub25', 'sub26', 'sub27', 'sub28', 'sub29', 'sub30', 'sub31', 'sub32'

% 定义所有条件
conditions = {'Need_nature', 'Noneed_nature', 'Freedom_nature', ...
             'Need_strategy', 'Noneed_strategy', 'Freedom_strategy'};

% 读取听觉onsets
onsets_auditory = struct();
for i = 1:length(auditory_subjects)
    subject = auditory_subjects{i};
    
    % 读取每个sub对应工作表的数据
    filepath = 'F:\2025新实验预处理\onsets.xlsx';
    sheet = subject;
    data = readtable(filepath, 'Sheet', sheet);
    
    % 读取6个条件的onsets（假设每个条件对应Excel中的一列）
    for c = 1:length(conditions)
        condition = conditions{c};
        condition_onsets = data{:, c};
        
        % 去除 NaN 值
        condition_onsets = condition_onsets(~isnan(condition_onsets));
        
        % 保存到结构中
        onsets_auditory.(subject).(condition) = condition_onsets;
        
        % 设置固定的持续时间
        durations_auditory.(subject).(condition) = repmat(fixed_duration, size(condition_onsets));
    end
end

% 保存onsets和durations为.mat文件
save('F:\2025新实验预处理\onsets_auditory.mat', 'onsets_auditory', 'durations_auditory');
