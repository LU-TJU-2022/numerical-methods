fid = fopen('F:\����\��\��������\case3\model4ʱ�䲽12000���10000����case3.txt');
a = fscanf(fid, '%f',[11 50005]);
fclose(fid);
x6=a(1,:); %��ȡ��һ�����ݣ�����x����
y6=a(2,:); %��ȡ�ڶ������ݣ�����y����
z6=a(3,:); %��ȡ���������ݣ����ӵ�Ũ��

[X,Y,Z]=griddata(x6,y6,z6,linspace(-6400, 6400,129)',linspace(-6400,6400,129),'linear');%��ֵ
pcolor(X,Y,Z);shading interp%α��ɫͼ
figure(10)
axis([474000 490000 3484000 3510000])

scatter(a(1,:),a(2,:),10,z6,'filled');
figure(20)
mesh(X,Y,Z)
contourf(X,Y,Z,10)

%---------y=3200������Ա�--------%
figure(30)
plot(X(1,:),Z(97,:), '-o')%SPH�����ֵ����
hold on
plot(X(1,:),c(:,101))
title('dx=0.01 dt=0.01')