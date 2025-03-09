fid = fopen('F:\研三\书\综述算例\case1\SPH时间步560结果case1.txt');
a = fscanf(fid, '%f',[5 50005]);
fclose(fid);
x6=a(1,:); %提取第一行数据，粒子x坐标
y6=a(2,:); %提取第二行数据，粒子y坐标
z6=a(3,:); %提取第三行数据，粒子的浓度

[X,Y,Z]=griddata(x6,y6,z6,linspace(-3,3,65)',linspace(-3,3,65),'linear');%插值
pcolor(X,Y,Z);shading interp%伪彩色图
figure(10)
axis([474000 490000 3484000 3510000])

scatter(a(1,:),a(2,:),10,z6,'filled');
figure(20)
mesh(X,Y,Z)
contourf(X,Y,Z,2)
%contourf(X,Y,Z,13)
%b=load('边界粒子.csv');
%hold on
%scatter(b(:,1),b(:,2),10,'filled','r')
%plot(X(1,:),Z(101,:))

%---------顶点处截面对比--------%
figure(30)
plot(X(1,:),Z(121,:), '-o')%SPH结果峰值截面
hold on
plot(X(1,:),c(:,101))
title('dx=0.01 dt=0.01')