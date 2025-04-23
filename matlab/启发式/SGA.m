%% 标准遗传算法SGA
clear
clc
pc = 0.7;
pm = 0.05;
% 定义遗传算法参数
NIND = 40;                                          % 个体数目
MAXGEN = 5000;                                       % 最大遗传次数
NVAR = 2;                                           % 变量的维数
PRECI = 20;                                         % 变量的二进制位数
GGAP = 0.9;                                         % 代沟
trace = zeros(MAXGEN, 1);
FieldD = [repmat(PRECI,[1,NVAR]); [-3,4.1; 12.1,5.8]; repmat([1;0;1;1],[1,NVAR])];    % 区域描述器
Chrom = crtbp(NIND, NVAR*PRECI);                    % 创建初始种群
gen = 0;                                            % 代计数器
maxY = 0;                                           % 最优值
ObjV = ObjectFunction(bs2rv(Chrom, FieldD));        % 计算初始种群个体的目标函数值
while gen < MAXGEN                                  % 迭代
    FitnV = ranking(-ObjV);                         % 分配适应度值
    SelCh = select('sus', Chrom, FitnV, GGAP);      % 选择
    SelCh = recombin('xovsp', SelCh, pc);           % 重组
    SelCh = mut(SelCh, pm);                         % 变异
    ObjVSel = ObjectFunction(bs2rv(SelCh, FieldD)); % 计算子代目标函数值
    [Chrom ObjV] = reins(Chrom, SelCh, 1, 1, ObjV, ObjVSel);    % 重插入
    gen = gen + 1;
    if maxY < max(ObjV)
        maxY = max(ObjV);
    end
    trace(gen, 1) = maxY;
end
%% 进化过程图
plot(1: gen, trace(:, 1));
%% 输出最优解
[Y, I] = max(ObjV);
X = bs2rv(Chrom, FieldD);
disp(['最优值为：', num2str(Y)]);
disp(['对应的自变量取值：', num2str(X(I, :))]);