fid = fopen('F:\����\��\��������\case1\SPHʱ�䲽560���case1.txt');
a = fscanf(fid, '%f',[5 50005]);
fclose(fid);
x6=a(1,:); %��ȡ��һ�����ݣ�����x����
y6=a(2,:); %��ȡ�ڶ������ݣ�����y����
z6=a(3,:); %��ȡ���������ݣ����ӵ�Ũ��

[X,Y,Z]=griddata(x6,y6,z6,linspace(-3,3,65)',linspace(-3,3,65),'linear');%��ֵ
pcolor(X,Y,Z);shading interp%α��ɫͼ
figure(10)
axis([474000 490000 3484000 3510000])

scatter(a(1,:),a(2,:),10,z6,'filled');
figure(20)
mesh(X,Y,Z)
contourf(X,Y,Z,2)
%contourf(X,Y,Z,13)
%b=load('�߽�����.csv');
%hold on
%scatter(b(:,1),b(:,2),10,'filled','r')
%plot(X(1,:),Z(101,:))

%---------���㴦����Ա�--------%
figure(30)
plot(X(1,:),Z(121,:), '-o')%SPH�����ֵ����
hold on
plot(X(1,:),c(:,101))
title('dx=0.01 dt=0.01')