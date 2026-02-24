% 定义NIfTI文件的文件夹路径
fmri_dir = 'E:\PFC\TMS\sub26\preproc\sub26';
output_filename = 'combined_4d_data.nii';  % 最终输出的4D NIfTI文件名

% 获取文件夹中的所有NIfTI文件名
file_list = dir(fullfile(fmri_dir, '*.nii'));
if isempty(file_list)
    file_list = dir(fullfile(fmri_dir, '*.nii.gz'));
end

% 如果没有找到NIfTI文件，报错并退出
if isempty(file_list)
    error('没有找到NIfTI文件！');
end

% 读取第一个NIfTI文件，获取图像尺寸
first_img = load_untouch_nii(fullfile(fmri_dir, file_list(1).name));
[x, y, z] = size(first_img.img);  % 获取图像的三维尺寸

n_files = length(file_list);
batch_size = 50;  % 更小的批次大小
num_batches = ceil(n_files/batch_size);

% 创建输出NIfTI文件的模板
template_nii = first_img;
template_nii.hdr.dime.dim(1) = 4;  % 设置为4D
template_nii.hdr.dime.dim(5) = n_files;  % 设置时间点数量

% 创建最终文件的占位符（一个空的3D volume）
template_nii.img = zeros(x, y, z, 1);
save_untouch_nii(template_nii, fullfile(fmri_dir, output_filename));

% 打开文件以进行二进制写入
fid = fopen(fullfile(fmri_dir, output_filename), 'r+');
if fid == -1
    error('无法打开输出文件进行写入');
end

% 获取数据开始的偏移量（跳过头文件）
initial_offset = template_nii.hdr.dime.vox_offset;
fseek(fid, initial_offset, 'bof');

% 按批次处理并直接写入文件
for batch = 1:num_batches
    start_idx = (batch-1)*batch_size + 1;
    end_idx = min(batch*batch_size, n_files);
    curr_batch_size = end_idx - start_idx + 1;
    
    % 为当前批次创建较小的临时数组
    batch_data = zeros(x, y, z, curr_batch_size);
    
    % 读取当前批次的数据
    for i = 1:curr_batch_size
        file_idx = start_idx + i - 1;
        current_img = load_untouch_nii(fullfile(fmri_dir, file_list(file_idx).name));
        batch_data(:,:,:,i) = current_img.img;
    end
    
    % 计算在文件中的位置
    offset = initial_offset + (start_idx-1) * x * y * z * 4; % 假设单精度浮点数（4字节）
    fseek(fid, offset, 'bof');
    
    % 将数据写入文件
    fwrite(fid, batch_data(:), 'float32');
    
    % 清理内存
    clear batch_data current_img
    
    fprintf('已处理 %d/%d 个文件\n', end_idx, n_files);
end

% 关闭文件
fclose(fid);

disp('4D数据已保存。');