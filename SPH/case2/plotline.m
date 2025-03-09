fid = fopen('F:\研三\书\综述算例\case2\model2时间步3500结果case2.txt');
a = fscanf(fid, '%f',[5 50005]);
fclose(fid);
x6=a(1,:); %提取第一行数据，粒子x坐标
y6=a(2,:); %提取第二行数据，粒子y坐标
z6=a(3,:); %提取第三行数据，粒子的浓度

[X,Y,Z]=griddata(x6,y6,z6,linspace(0, 6400,129)',linspace(0,6400,129),'linear');%插值
pcolor(X,Y,Z);shading interp%伪彩色图
figure(10)
axis([474000 490000 3484000 3510000])

scatter(a(1,:),a(2,:),10,z6,'filled');
figure(20)
mesh(X,Y,Z)
contourf(X,Y,Z,10)

%---------y=3200处截面对比--------%
figure(30)
plot(X(1,:),Z(97,:), '-o')%SPH结果峰值截面
hold on
plot(X(1,:),c(:,101))
title('dx=0.01 dt=0.01')