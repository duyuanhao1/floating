function Objective = powercalculate(wind,address,site)
%powercalculate 计算风电场功率，输入为风机位置与风向、风速，
% 格式为wind(v1^3 ang1 p1;v2^3 ang2 p2;.....)，address(x1 y1;x2 y2...）
% 输出为风电场功率power
% load semi.mat
% load TLP.mat
% load spar.mat

% 为每种风机平台定义不同的概率矩阵
prob_matrix_y1=[0.5/8,0.5/8,0.5/8,0.5/8,0.5,0.5/8,0.5/8,0.5/8,0.5/8]';
% prob_matrix_y1 = prob_matrix_y1 / sum(prob_matrix_y1); % 归一化
prob_matrix_y2=[0,0,0,0,1,0,0,0,0]';
% prob_matrix_y2 = prob_matrix_y2 / sum(prob_matrix_y2); % 归一化
prob_matrix_y3=[0.2/8,0.2/8,0.2/8,0.2/8,0.2,0.2/8,0.2/8,0.2/8,0.2/8]';
% prob_matrix_y3 = prob_matrix_y3 / sum(prob_matrix_y3); % 归一化

% 调用函数计算尾流影响
DD1_1 = calculate_wake_effect_3(address, wind, prob_matrix_y1,prob_matrix_y1);%y1风机平台对y2风机平台的尾流影响
DD1_2 = calculate_wake_effect_3(address, wind, prob_matrix_y1,prob_matrix_y2);
DD1_3 = calculate_wake_effect_3(address, wind, prob_matrix_y1,prob_matrix_y3);

DD2_1 = calculate_wake_effect_3(address, wind, prob_matrix_y2,prob_matrix_y1);
DD2_2 = calculate_wake_effect_3(address, wind, prob_matrix_y2,prob_matrix_y2);
DD2_3 = calculate_wake_effect_3(address, wind, prob_matrix_y2,prob_matrix_y3);

DD3_1 = calculate_wake_effect_3(address, wind, prob_matrix_y3,prob_matrix_y1);
DD3_2 = calculate_wake_effect_3(address, wind, prob_matrix_y3,prob_matrix_y2);
DD3_3 = calculate_wake_effect_3(address, wind, prob_matrix_y3,prob_matrix_y3);
% load DD11_1.mat;
% load DD12_1.mat;
% load DD13_1.mat;
% load DD21_1.mat;
% load DD22_1.mat;
% load DD23_1.mat;
% load DD31_1.mat;
% load DD32_1.mat;
% load DD33_1.mat;
I = size(address, 1); % I为将要选取的地址的离散集合,格式为address(x1 y1;x2 y2...）
D = size(wind, 1); % D为风的离散集合,格式为wind(v1 ang1 p1;v2 ang2 p2;.....），角度为与正北方向按顺时针计算的角度大小

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

% DD=Wake_effect_v3(address,wind);
% WW=DD.^2;  %vv为对应的速度比例
I=size(address,1); %I为将要选取的地址的离散集合,格式为address(x1 y1;x2 y2...）
D=size(wind,1);%D为风的离散集合,格式为wind(v1 ang1 p1;v2 ang2 p2;.....），角度为与正北方向按顺时针计算的角度大小
power_rate=zeros(size(address,1),size(wind,1));
power_temp=0;
M=5000;diameter=40;constraintviolation=0;
for j=1:I
       for k=1:j-1
          if(norm(address(j,:)-address(k,:))<5*diameter) %安全距离为5D    
              constraintviolation=constraintviolation+1;
              break %跳出第一层循环
          end
       end
%       if(constraintviolation~=0) %安全距离为5D    
%           break  %跳出第二层循环
%       end
end

%% 功率验算
y1=site(:,1);
y2=site(:,2);
y3=site(:,3);
site1=sum(site,2);
loss = 0;
power = 0;
for dd = 1:D
    for ii = 1:I
        for jj = 1:I           
            loss = loss + ((y1(ii)) * ((y1(jj)) * loss_num1_1(ii, jj, dd) + (y2(jj)) * loss_num1_2(ii, jj, dd) + (y3(jj)) * loss_num1_3(ii, jj, dd)) + ...
                         (y2(ii)) * ((y1(jj)) * loss_num2_1(ii, jj, dd) + (y2(jj)) * loss_num2_2(ii, jj, dd) + (y3(jj)) * loss_num2_3(ii, jj, dd)) + ...
                         (y3(ii)) * ((y1(jj)) * loss_num3_1(ii, jj, dd) + (y2(jj)) * loss_num3_2(ii, jj, dd) + (y3(jj)) * loss_num3_3(ii, jj, dd)));
        end
        power = power + 2.8935 * wind(dd, 1) * wind(dd, 3) * (1 - loss) * site1(ii) - M*constraintviolation;
        % fprintf('i 的值为 %d，j 的值为 %d,d为%d\n', ii, jj,dd);
        loss = 0;
    end 
end
p_rate=5000;
wake_loss = 3000*0.93*25*(p_rate*(sum(y1)+sum(y2)+sum(y3))-power);
%wake_loss=-sum(p,'all');
fp_cost = sum(y1)*246800000+sum(y2)*288100000+sum(y3)*282200000;
Objective = wake_loss + fp_cost;
end