function FitnV = ranking(ObjV)
    % ObjV: 目标函数值矩阵，每一行代表一个个体的目标函数值
    % FitnV: 适应度值向量，每个元素对应一个个体的适应度
    
    [Nind, ~] = size(ObjV);
    
    % 对目标函数值进行排序，得到排序后的索引
    [~, sortedIndices] = sort(ObjV);
    
    % 线性排序法，适应度值从 1 到 Nind 线性递增
    FitnV = zeros(Nind, 1);
    for i = 1:Nind
        FitnV(sortedIndices(i)) = i;
    end
    
    % 归一化适应度值，使其总和为 Nind
    FitnV = FitnV * Nind / sum(FitnV);
end