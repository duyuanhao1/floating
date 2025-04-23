function obj = ObjectFunction(X)
%% 待优化的目标函数
% X 的每行为一个个体,每行格式为[x1 x2 x3 x4……xn y1 y2…… yn]
col = size(X, 1);
wind=[];
for n = 0.5:1:11.5
    wind = [wind; [12^3, 2*pi/12*(n-0.5)+pi, 1/12]];
end
% wind=[12^3,0,1];
% load('winddata_case3.mat', 'wind');
for i = 1: col
    address=reshape(X(i,:),[],2);
    site= readmatrix('4_12_101024.txt');  % 赋初值
    obj(i, 1) = powercalculate(wind,address,site);
end
end