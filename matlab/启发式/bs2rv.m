function phen = bs2rv(Chrom, FieldD)
    % 此函数用于将二进制编码转换为实数值
    % Chrom: 二进制编码的种群
    % FieldD: 译码矩阵
    
    [Nind, Lind] = size(Chrom);
    [~, NVars] = size(FieldD);
    
    BaseV = 2.^([FieldD(1, 1)-1:-1:0]);
    BaseV = repmat(BaseV, Nind, 1);
    
    % 提取每个变量的二进制位数
    Prec = FieldD(1, :);
    % 提取每个变量的上下界
    LB = FieldD(2, :);
    UB = FieldD(3, :);
    
    phen = zeros(Nind, NVars);
    
    start = 1;
    for i = 1:NVars
        % 提取当前变量的二进制编码
        subChrom = Chrom(:, start:start+Prec(i)-1);
        % 将二进制编码转换为十进制数
        decVal = sum(subChrom.*BaseV(:, 1:Prec(i)), 2);
        % 线性映射到实际取值范围
        phen(:, i) = LB(i) + decVal * (UB(i) - LB(i)) / (2^Prec(i) - 1);
        start = start + Prec(i);
    end
end