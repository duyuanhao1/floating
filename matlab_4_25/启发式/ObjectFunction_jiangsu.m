function obj = ObjectFunction_jiangsu(X)
%% 待优化的目标函数
% X 的每行为一个个体,每行格式为[x1 x2 x3 x4……xn y1 y2…… yn]
col = size(X, 1);
load ('wind_data.mat','prob');
%prob 为风速对应的概率，分别为12m/s及以上,11m/s(到12m/s),10m/s,9m/s,8,7,6,5,4,3,0-3
v_tri=[11.5:-1:3.5].^3';
wind=[];
for n=1:36
    wind=[wind;[prob(n,2:10)*v_tri,2*pi/36*(n-10),sum(prob(n,2:10))]];
end 
for i = 1: col
    address=reshape(X(i,:),[],2);
    obj(i, 1) = powercalculate_jiangsu(wind,address,prob);
end
end