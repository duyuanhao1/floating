%% ����Ⱥ�Ŵ��㷨
clear;
global NVAR % ����ȫ�ֱ��� NVAR
NT=40;
NIND = 40;                                  % ������Ŀ
NVAR = 2*NT;                                % ������ά��
PRECI = 20;                                 % �����Ķ�����λ��
GGAP = 0.9;                                 % ����
MP = 4;                                    % ��Ⱥ��Ŀ
start1=tic;
address=[];
gridlength=580;gridwidth=580;
for x=0:gridlength:5220
    for y=0:gridwidth:5220
         address=[address;[x y]];
    end
end
load wind_data.mat
%prob Ϊ���ٶ�Ӧ�ĸ��ʣ��ֱ�Ϊ12m/s������,11m/s(��12m/s),10m/s,9m/s,8,7,6,5,4,3,0-3
v_tri=[11.5:-1:3.5].^3';
wind=[];
for n=1:36
    wind=[wind;[prob(n,2:10)*v_tri,2*pi/36*(n-10),sum(prob(n,2:10))]];
end 
site=readmatrix('quarcombine_for_jiangsu_owf.txt');
bound = sectiongeneration(address,site,gridlength,gridwidth);
disp(['Size of bound: ', num2str(size(bound))]); % ��ʾ bound �Ĵ�С
disp(['NVAR: ', num2str(NVAR)]); % ��ʾ NVAR ��ֵ
FieldD = [repmat(PRECI,[1,NVAR]);bound; repmat([1;0;1;1],[1,NVAR])]; %  �������
for i=1:MP
    Chrom{i} = crtbp(NIND, NVAR*PRECI);     % ������ʼ��Ⱥ
end
pc = 0.7 + (0.9 - 0.7)*rand(MP, 1);         % ��[0.7,0.9]��Χ����������������
pm = 0.001 + (0.05 - 0.001)*rand(MP, 1);    % ��[0.001,0.05]��Χ����������������
gen = 0;                                    % ��ʼ�Ŵ�����
gen0 = 0;                                   % ��ʼ���ִ���
MAXGEN = 100;                              % ���Ÿ������ٱ��ִ���
maxY = 0;                                   % ����ֵ
for i = 1: MP
    ObjV{i} = ObjectFunction_jiangsu(bs2rv(Chrom{i}, FieldD));                      % �������ʼ��Ⱥ�����Ŀ�꺯��ֵ
end
MaxObjV = zeros(MP, 1);                     % ��¼������Ⱥ
MaxChrom = zeros(MP, PRECI*NVAR);           % ��¼������Ⱥ�ı���
while gen0 <= MAXGEN
    gen = gen + 1;                          % �Ŵ�������1
    for i = 1: MP
        FitnV{i} = ranking(-ObjV{i});                           % ����Ⱥ����Ӧ��
        SelCh{i} = select('sus', Chrom{i}, FitnV{i}, GGAP);     % ѡ�����
        SelCh{i} = recombin('xovsp', SelCh{i}, pc(i));          % �������
        SelCh{i} = mut(SelCh{i}, pm(i));                        % �������
        ObjVSel = ObjectFunction_jiangsu(bs2rv(SelCh{i}, FieldD));      % �����Ӵ�Ŀ�꺯��ֵ
        [Chrom{i}, ObjV{i}] = reins(Chrom{i}, SelCh{i}, 1, 1, ObjV{i}, ObjVSel);    % �ز������
    end
    [Chrom, ObjV] = immigrant(Chrom,ObjV);  % �������
    [MaxObjV, MaxChrom] = EliteInduvidual(Chrom, ObjV, MaxObjV, MaxChrom);          % �˹�ѡ�񾫻���Ⱥ
    YY(gen) = max(MaxObjV);                 % �ҳ�������Ⱥ�����ŵĸ���
    if YY(gen) > maxY+10                       % �жϵ�ǰ�Ż�ֵ�Ƿ���ǰһ���Ż�ֵ��ͬ
        maxY = YY(gen);                     % ��������ֵ
        gen0 = 0;
    else
        gen0 = gen0+1;                      % ����ֵ���ִ�����1
    end
    
end
time1=toc(start1)
%% ��������ͼ
figure(1)
plot(1: gen, max(0,YY))
xlabel('��������')
ylabel('���Ž�仯')
title('��������')
xlim([1, gen])
%% ������Ž�
[Y,I] = max(MaxObjV);    % �ҳ�������Ⱥ�����ŵĸ���
X = bs2rv(MaxChrom(I, :), FieldD);   % ���Ÿ���Ľ����
disp(['����ֵΪ��', num2str(Y)])
disp(['��Ӧ���Ա���ȡֵ��', num2str(X)])
realad=reshape(X,[],2);
figure(2)
for i=1:size(realad,1)
    plot(realad(i,1),realad(i,2),'vk','LineWidth',2,'MarkerFaceColor','k');
    hold on
end
ax=gca;
% ax.YTick = [0:580:5800];
% ax.XTick = [0:580:5800];
% ax.XLim = [0,5800];
% ax.YLim = [0,5800];
grid on
grid on
set(gcf,'color','white');%��ɫ
set(gca,'LineWidth',2);
set(gca,'GridAlpha',1);
title('���΢��ѡַ���ַ���');xlabel('����λ��/m');ylabel('����λ��/m');
%% �����糡����
power=powercalculate_jiangsu(wind,realad,prob)