%% 多种群遗传算法
clear
close all
NT=40;
NIND = 30;                                  % 个体数目
NVAR = 2*NT;                                % 变量的维数
PRECI = 20;                                 % 变量的二进制位数
GGAP = 0.9;                                 % 代沟
MP = 5;                                    % 种群数目
wind=[];
for n = 0.5:1:35.5
    wind = [wind; [12^3, 2*pi/36*(n-0.5)+pi, 1/36]];
end
address=[];gridlength=200;gridwidth=200;
for x=0:gridlength:1800
    for y=0:gridwidth:1800
         address=[address;[x y]];
    end
end
site=readmatrix('4_24_9101036.txt');
site=sum(site,2);
bound = sectiongeneration_v2(address,site,gridlength,gridwidth);
FieldD = [repmat(PRECI,[1,NVAR]);bound; repmat([1;0;1;1],[1,NVAR])]; %  译码矩阵
for i=1:MP
    Chrom{i} = crtbp(NIND, NVAR*PRECI);     % 创建初始种群
end
pc = 0.7 + (0.9 - 0.7)*rand(MP, 1);         % 在[0.7,0.9]范围内随机产生交叉概率
pm = 0.001 + (0.05 - 0.001)*rand(MP, 1);    % 在[0.001,0.05]范围内随机产生变异概率
gen = 0;                                    % 初始遗传代数
gen0 = 0;                                   % 初始保持代数
MAXGEN = 100;                                % 最优个体最少保持代数
minY = 0;                                   % 最优值
for i = 1: MP
    ObjV{i} = ObjectFunction(bs2rv(Chrom{i}, FieldD));                      % `计算各初始种群个体的目标函数值
end
MaxObjV = zeros(MP, 1);                     % 记录精华种群
MaxChrom = zeros(MP, PRECI*NVAR);           % 记录精华种群的编码
while gen0 <= MAXGEN
    gen = gen + 1;                          % 遗传代数加1
    for i = 1: MP
        FitnV{i} = ranking(-ObjV{i});                           % 各种群的适应度
        SelCh{i} = select('sus', Chrom{i}, FitnV{i}, GGAP);     % 选择操作
        SelCh{i} = recombin('xovsp', SelCh{i}, pc(i));          % 交叉操作
        SelCh{i} = mut(SelCh{i}, pm(i));                        % 变异操作
        ObjVSel =  ObjectFunction(bs2rv(SelCh{i}, FieldD));      % 计算子代目标函数值
        [Chrom{i}, ObjV{i}] = reins(Chrom{i}, SelCh{i}, 1, 1, ObjV{i}, ObjVSel);    % 重插入操作
    end
    [Chrom, ObjV] = immigrant(Chrom,ObjV);  % 移民操作
    [MaxObjV, MaxChrom] = EliteInduvidual(Chrom, ObjV, MaxObjV, MaxChrom);          % 人工选择精华种群
    YY(gen) = min(MaxObjV);                 % 找出精华种群中最优的个体
    if YY(gen) < minY                       % 判断当前优化值是否与前一次优化值相同
        minY = YY(gen);                     % 更新最优值
        gen0 = 0;
    else
        gen0 = gen0+1;                      % 最优值保持次数加1
    end
end
%% 进化过程图
figure(1)
plot(1: gen, max(0,YY),'LineWidth',2)
xlabel('evolutional generation')
ylabel('fitness')
% title('进化过程')
xlim([1, gen])
set(gcf,'color','white');%白色
set(gca,'LineWidth',1);
%% 输出最优解
[Y,I] = max(MaxObjV);    % 找出精华种群中最优的个体
X = (bs2rv(MaxChrom(I, :), FieldD));   % 最优个体的解码解
disp(['最优值为：', num2str(Y)])
disp(['对应的自变量取值：', num2str(X)])
realad=reshape(X,[],2);
figure(2)
for i=1:size(realad,1)
    plot(realad(i,1),realad(i,2),'vk','LineWidth',4,'MarkerFaceColor','k');
    hold on
end
grid on
set(gcf,'color','white');%白色
set(gca,'LineWidth',2);
set(gca,'GridAlpha',1);
% title('风机微观选址布局方案');xlabel('横轴位置/m');ylabel('纵轴位置/m');
%% 验算风电场功率

money=powercalculate(wind,realad,site)