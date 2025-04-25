clear;yalmip('clear');
% 初始化地址和风向矩阵
address = [];
wind = [];
gn = 9; % (9)+1×(9)+1的网格
gridlength = 200;
gridwidth = 200;

% 生成地址矩阵
for x = 0:gridlength:(gridwidth*gn)
    for y = 0:gridwidth:(gridwidth*gn)
         address = [address; [x y]];
    end
end
for n = 0.5:1:35.5
    wind = [wind; [12^3, 2*pi/36*(n-0.5)+pi, 1/36]];
end
% wind = [12^3, 0, 1];
p_rate = 5000;

% 为每种风机平台定义不同的概率矩阵
prob_matrix_y1=[0.5/8,0.5/8,0.5/8,0.5/8,0.5,0.5/8,0.5/8,0.5/8,0.5/8]';
% prob_matrix_y1 = prob_matrix_y1 / sum(prob_matrix_y1); % 归一化
prob_matrix_y2=[0,0,0,0,1,0,0,0,0]';
% prob_matrix_y2 = prob_matrix_y2 / sum(prob_matrix_y2); % 归一化
prob_matrix_y3=[0.8/8,0.8/8,0.8/8,0.8/8,0.2,0.8/8,0.8/8,0.8/8,0.8/8]';
% prob_matrix_y3 = prob_matrix_y3 / sum(prob_matrix_y3); % 归一化
 load DD1_1_36_12.mat;
 load DD1_2_36_12.mat;
 load DD1_3_36_12.mat;
 load DD2_1_36_12.mat;
 load DD2_2_36_12.mat;
 load DD2_3_36_12.mat;
 load DD3_1_36_12.mat;
 load DD3_2_36_12.mat;
 load DD3_3_36_12.mat;

% DD1_1 = calculate_wake_effect_3(address, wind, prob_matrix_y1,prob_matrix_y1);%y1风机平台对y2风机平台的尾流影响
% DD1_2 = calculate_wake_effect_3(address, wind, prob_matrix_y1,prob_matrix_y2);
% DD1_3 = calculate_wake_effect_3(address, wind, prob_matrix_y1,prob_matrix_y3);
% 
% DD2_1 = calculate_wake_effect_3(address, wind, prob_matrix_y2,prob_matrix_y1);
% DD2_2 = calculate_wake_effect_3(address, wind, prob_matrix_y2,prob_matrix_y2);
% DD2_3 = calculate_wake_effect_3(address, wind, prob_matrix_y2,prob_matrix_y3);
% 
% DD3_1 = calculate_wake_effect_3(address, wind, prob_matrix_y3,prob_matrix_y1);
% DD3_2 = calculate_wake_effect_3(address, wind, prob_matrix_y3,prob_matrix_y2);
% DD3_3 = calculate_wake_effect_3(address, wind, prob_matrix_y3,prob_matrix_y3);

I = size(address, 1);
D = size(wind, 1); 

vv1_1 = 1 - DD1_1;  % vv为对应的速度比例
vv1_2 = 1 - DD1_2;
vv1_3 = 1 - DD1_3;
vv2_1 = 1 - DD2_1;
vv2_2 = 1 - DD2_2;
vv2_3 = 1 - DD2_3;
vv3_1 = 1 - DD3_1;
vv3_2 = 1 - DD3_2;
vv3_3 = 1 - DD3_3;

loss_num1_1 = 1 - vv1_1.^3;
loss_num1_2 = 1 - vv1_2.^3;
loss_num1_3 = 1 - vv1_3.^3;
loss_num2_1 = 1 - vv2_1.^3;
loss_num2_2 = 1 - vv2_2.^3;
loss_num2_3 = 1 - vv2_3.^3;
loss_num3_1 = 1 - vv3_1.^3;
loss_num3_2 = 1 - vv3_2.^3;
loss_num3_3 = 1 - vv3_3.^3;

y11 = binvar(size(address, 1), 1, 'full'); % 半潜  
y22 = binvar(size(address, 1), 1, 'full'); % 张力腿 
y33 = binvar(size(address, 1), 1, 'full'); % 立柱 
y = binvar(size(address, 1), 1, 'full'); % 立柱 
p = sdpvar(size(address, 1), size(wind, 1), 'full');

%% 添加约束
Constraints = [sum(y11 + y22 + y33, "all") == 40;p >= 0];  % 约束，总机数

Constraints=[Constraints,y == y11 + y22 + y33];

for dd = 1:D
    Constraints = [Constraints, p(:, dd) <= 2.8935 * wind(dd, 3) * wind(dd, 1) * y]; % 最大风速约束
end

for dd = 1:D
    for i = 1:I
        prev_loss1 = 0;
        prev_loss2 = 0;
        prev_loss3 = 0;
        for j = 1:I
            prev_loss1 = prev_loss1 + y11(j) * loss_num1_1(i, j, dd) + y22(j) * loss_num1_2(i, j, dd) + y33(j) * loss_num1_3(i, j, dd);
            prev_loss2 = prev_loss2 + y11(j) * loss_num2_1(i, j, dd) + y22(j) * loss_num2_2(i, j, dd) + y33(j) * loss_num2_3(i, j, dd);
            prev_loss3 = prev_loss3 + y11(j) * loss_num3_1(i, j, dd) + y22(j) * loss_num3_2(i, j, dd) + y33(j) * loss_num3_3(i, j, dd);
            
        end

        Constraints = [Constraints,... 
            p(i, dd) <= 2.8935 * wind(dd, 3) * wind(dd, 1) * (y11(i)+y22(i)+y33(i)-prev_loss1*y11(i)-prev_loss2*y22(i)-prev_loss3*y33(i))];      
       
    end
end
 fprintf('i 的值为 %d，j 的值为 %d,d为%d\n', i, j ,dd);
% 添加新的距离约束,y1半潜40，y2张力腿0，y3立柱80
for i = 1:I%19800个
    for j = i + 1:I
        dist = norm(address(i, :) - address(j, :));
        % y1处的风机距离内不许有y1风机
        Constraints = [Constraints, y11(i) + y11(j) <= 1 + (dist >= 280)];
        % y2处风机距离内不许有y1风机
        Constraints = [Constraints, y11(i) + y22(j) <= 1 + (dist >= 240)];
        % y3处风机距离内不许有y3风机
        Constraints = [Constraints, y33(i) + y33(j) <= 1 + (dist >= 360)];
        % y3处风机距离内不许有y2风机
        Constraints = [Constraints, y33(i) + y22(j) <= 1 + (dist >= 280)];
        % y3处风机距离内不许有y1风机
        Constraints = [Constraints, y33(i) + y11(j) <= 1 + (dist >= 320)];
    end
end
wake_loss = 7300*0.93*30*(p_rate*(sum(y11)+sum(y22)+sum(y33))-sum(p,'all'));

fp_cost = sum(y11)*95600000+sum(y22)*111000000+sum(y33)*116900000;
Objective = wake_loss + fp_cost;

ops = sdpsettings('solver', 'gurobi','verbose', 2 ,'gurobi.NonConvex',2,'usex0',1);%设置求解器



%% 求解并输出结果
solution = optimize(Constraints, Objective, ops);

%% 绘制图像
site = [];
site(:, 1) = value(y11);
site(:, 2) = value(y22);
site(:, 3) = value(y33);
turbine_address = [];
turbine_address = [];
hold on;
for ii = 1:(gn+1)*(gn+1)
        if((0.999<=site(ii, 1)&&site(ii, 1)<=1.001))
         turbine_address = [turbine_address; address(ii, 1), address(ii, 2)];
         plot(address(ii, 1)+gridlength/2, address(ii, 2)+gridwidth/2, 'vk', 'LineWidth', 5, 'MarkerFaceColor', 'k'); % 三角半潜
         hold on;
        end
end
for ii = 1:(gn+1)*(gn+1)
        if((0.999<=site(ii, 2)&&site(ii, 2)<=1.001))
         turbine_address = [turbine_address; address(ii, 1), address(ii, 2)];
         plot(address(ii, 1)+gridlength/2, address(ii, 2)+gridwidth/2, 'ok', 'LineWidth', 6, 'MarkerFaceColor', 'w'); % 圆圈张力腿
         hold on;
        end
end
for ii = 1:(gn+1)*(gn+1)
        if((0.999<=site(ii, 3)&&site(ii, 3)<=1.001))
         turbine_address = [turbine_address; address(ii, 1), address(ii, 2)];
         plot(address(ii, 1)+gridlength/2, address(ii, 2)+gridwidth/2, 'sk', 'LineWidth', 6, 'MarkerFaceColor', 'k'); % 正方形星立柱
         hold on;
        end
end

ax = gca;
ax.YTick = [0:gridlength:gridlength*(gn+1)];
ax.XTick = [0:gridlength:gridlength*(gn+1)];
ax.XLim = [0, gridlength*(gn+1)];
ax.YLim = [0, gridlength*(gn+1)];
grid on;
% title('风机微观选址布局方案');xlabel('横轴位置/m');ylabel('纵轴位置/m');
set(gcf, 'color', 'white'); % 白色
set(gca, 'LineWidth', 2);
set(gca, 'GridAlpha', 1);
% title('风机微观选址布局方案');xlabel('横轴位置/m');ylabel('纵轴位置/m');

%% 功率验算1
% site1=readmatrix('4_24_9101036.txt');
% y11=site(:,1);
% y22=site(:,2);
% y33=site(:,3);
site1 = value(y);
y11=value(y11);
y22=value(y22);
y33=value(y33);
loss = 0;
power = 0;
for dd = 1:D
    for ii = 1:I
        for jj = 1:I           
            loss = loss + ((y11(ii)) * ((y11(jj)) * loss_num1_1(ii, jj, dd) + (y22(jj)) * loss_num1_2(ii, jj, dd) + (y33(jj)) * loss_num1_3(ii, jj, dd)) + ...
                         (y22(ii)) * ((y11(jj)) * loss_num2_1(ii, jj, dd) + (y22(jj)) * loss_num2_2(ii, jj, dd) + (y33(jj)) * loss_num2_3(ii, jj, dd)) + ...
                         (y33(ii)) * ((y11(jj)) * loss_num3_1(ii, jj, dd) + (y22(jj)) * loss_num3_2(ii, jj, dd) + (y33(jj)) * loss_num3_3(ii, jj, dd)));
        end
        power = power + 2.8935 * wind(dd, 1) * wind(dd, 3) * (1 - loss) * site1(ii);
        % fprintf('i 的值为 %d，j 的值为 %d,d为%d\n', ii, jj,dd);
        loss = 0;
    end 
end
power
power_obj = value(sum(p,'all'))
p_real = value(p);
v_real = (p_real/2.8935*36).^(1/3);
all_wake_loss = 7300*0.93*30*(p_rate*(sum(y11)+sum(y22)+sum(y33))-power);
%wake_loss=-sum(p,'all');
all_fp_cost = sum(y11)*95600000+sum(y22)*131000000+sum(y33)*136900000;
all_Objective = all_wake_loss + all_fp_cost;
Objective=value(Objective);