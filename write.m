function write(new_data)
    % 定义文件名
    filename = 'geo_example.txt';
    %读取原文件前三行
    fid = fopen(filename, 'r');
    header = cell(3, 1);
    for i = 1:3
        header{i} = fgetl(fid); % 逐行读取固定头信息
    end
    fclose(fid);
    
    %生成新数据数组
    % 假设新数据由函数生成，格式为 Nx2 的矩阵
   
    
    % 验证数据维度
    assert(size(new_data, 2) == 2, '新数据必须是Nx2的矩阵');
    assert(all(isnumeric(new_data(:))), '新数据必须为数值类型');
    
    %覆盖写入原文件
    fid = fopen(filename, 'w'); % 以写入模式清空文件
    
    % 写入头信息
    for i = 1:3
        fprintf(fid, '%s\n', header{i});
    end
    
    % 写入新数据（保持原始精度格式）
    formatSpec = '%.20e %.20e\n'; % 严格匹配原文件格式
    fprintf(fid, formatSpec, new_data'); % 注意转置矩阵
    
    fclose(fid);
    disp('新数据已成功写入原文件');
end