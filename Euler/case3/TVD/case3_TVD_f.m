function [] = case1_TVD_f(Limiter,method)
tic;
addpath('Limiters');
%% 数据初始�?
Xbegin=-6400;%//length of the domain
Xend=6400;
Ybegin=-6400;
Yend=6400;
T=12000;
Nx=128;
Ny=128;
M=12000;
dx=(Xend-Xbegin)/Nx;% //this is the size of the space step
dy=(Yend-Ybegin)/Ny;
dt=T/M; %//this is the size of the time step
D=1.0;%%%扩散系数
x1=zeros(1,Nx+1); % array allocation
y1=zeros(1,Ny+1);
c1=zeros(Nx+1,Ny+1);
u=zeros(Nx+1,Ny+1);
v=zeros(Nx+1,Ny+1);
x0=0;
y0=3200;
sigmax=707.15;
sigmay=500;

for i=1:Nx+1
    for j=1:Ny+1
    x1(i)=dx*(i-(Nx/2+1))+Xbegin; %//there is no cell 0 in scilab that is why we begin with 1
    y1(j)=dy*(j-(Ny/2+1))+Ybegin;
    u(i,j)=-pi/6000*y1(j);
    v(i,j)=pi/6000*x1(i);
    
    c1(i,j)=exp(-(x1(i)-x0).^2/2/sigmax/sigmax-(y1(j)-y0).^2/2/sigmay/sigmay);
    end
end

c2=c1;
[X,Y]=meshgrid(1:Nx+1,1:Ny+1);%%%%meshgrid,注意矩阵行数对应y，列数对应X，比�?*4列，则y=1:3，x=1:4

%% 主体计算流程
for n=1:M
    c1(1,:)=0;
    c1(Nx+1,:)=0;
    c1(:,Ny+1)=0;
    c1(:,1)=0;
    
    c1(2,:)=0;
    c1(Nx,:)=0;
    c1(:,Ny)=0;
    c1(:,2)=0;
    for i=3:Nx-1
        for j=3:Ny-1
        %%-------选择限制�?--------%%
        switch Limiter
            case 'Superbee'
                [Limiter1x, Limiter1y] = Superbee(c1, dx, dy, i, j);
            case 'Minmod'
                [Limiter1x, Limiter1y] = Minmod(c1, dx, dy, i, j);
            case 'Woodward'
                [Limiter1x, Limiter1y] = Woodward(c1, dx, dy, i, j);
            case 'Van Leer'
                [Limiter1x, Limiter1y] = VanLeer(c1, dx, dy, i, j);
            case 'UMIST'
                [Limiter1x, Limiter1y] = UMIST(c1, dx, dy, i, j);
            case 'WACEB'
                [Limiter1x, Limiter1y] = WACEB(c1, dx, dy, i, j);
            case 'Albada'
                [Limiter1x, Limiter1y] = Albada(c1, dx, dy, i, j);
            case 'OSPRE'
                [Limiter1x, Limiter1y] = OSPRE(c1, dx, dy, i, j);
            case 'TCDF'
                [Limiter1x, Limiter1y] = TCDF(c1, dx, dy, i, j);
            case 'Koren'
                [Limiter1x, Limiter1y] = Koren(c1, dx, dy, i, j);
        end

        Cx_i_min_1_2_R=c1(i,j)-0.5*dx*Limiter1x(1);%%c(i-1/2,R,x)
        Cy_i_min_1_2_R=c1(i,j)-0.5*dy*Limiter1y(1);%%c(i-1/2,R,x)
        Cx_i_plu_1_2_L=c1(i,j)+0.5*dx*Limiter1x(1);%%c(i-1/2,R,x)
        Cy_i_plu_1_2_L=c1(i,j)+0.5*dy*Limiter1y(1);%%c(i-1/2,R,x)
        Cx_i_min_1_2_L=c1(i-1,j)+0.5*dx*Limiter1x(2);%%c(i-1/2,R,x)
        Cy_i_min_1_2_L=c1(i,j-1)+0.5*dy*Limiter1y(2);%%c(i-1/2,R,x)
        Cx_i_plu_1_2_R=c1(i+1,j)-0.5*dx*Limiter1x(3);%%c(i-1/2,R,x)
        Cy_i_plu_1_2_R=c1(i,j+1)-0.5*dy*Limiter1y(3);%%c(i-1/2,R,x) 
        
        %%-------计算对流项F（i+1/2�?F（i-1/2�?--------%%
        %%-------计算f(Ci+1/2,R),f(Ci+1/2,L),F(i+1/2),可更换�?�?--------%%
        fc_x_i_plu_1_2_R=u(i,j)*Cx_i_plu_1_2_R;
        fc_y_i_plu_1_2_R=v(i,j)*Cy_i_plu_1_2_R;
        fc_x_i_plu_1_2_L=u(i,j)*Cx_i_plu_1_2_L;
        fc_y_i_plu_1_2_L=v(i,j)*Cy_i_plu_1_2_L;
        fc_x_i_min_1_2_R=u(i-1,j)*Cx_i_min_1_2_R;
        fc_y_i_min_1_2_R=v(i,j-1)*Cy_i_min_1_2_R;
        fc_x_i_min_1_2_L=u(i-1,j)*Cx_i_min_1_2_L;
        fc_y_i_min_1_2_L=v(i,j-1)*Cy_i_min_1_2_L;

        switch method
            case 'MUSCL'
                %%%%%%%%%%%%%%%%计算a（i+1/2�?
                if fc_x_i_plu_1_2_R ~= fc_x_i_plu_1_2_L
                    ax_i_plu_1_2_RL=(fc_x_i_plu_1_2_R-fc_x_i_plu_1_2_L)/(Cx_i_plu_1_2_R-Cx_i_plu_1_2_L);
                else
                    ax_i_plu_1_2_RL=u(i,j)*Cx_i_plu_1_2_R;
                end
                if fc_y_i_plu_1_2_R ~= fc_y_i_plu_1_2_L
                    ay_i_plu_1_2_RL=(fc_y_i_plu_1_2_R-fc_y_i_plu_1_2_L)/(Cy_i_plu_1_2_R-Cy_i_plu_1_2_L);
                else
                    ay_i_plu_1_2_RL=v(i,j)*Cy_i_plu_1_2_R;
                end
                %%%%%%%%%%%%%%%%计算a（i-1/2�?
                if fc_x_i_min_1_2_R ~= fc_x_i_min_1_2_L
                    ax_i_min_1_2_RL=(fc_x_i_min_1_2_R-fc_x_i_min_1_2_L)/(Cx_i_min_1_2_R-Cx_i_min_1_2_L);
                else
                    ax_i_min_1_2_RL=u(i,j)*Cx_i_min_1_2_R;
                end
                if fc_y_i_min_1_2_R ~= fc_y_i_min_1_2_L
                    ay_i_min_1_2_RL=(fc_y_i_min_1_2_R-fc_y_i_min_1_2_L)/(Cy_i_min_1_2_R-Cy_i_min_1_2_L);
                else
                    ay_i_min_1_2_RL=v(i,j)*Cy_i_min_1_2_R;
                end
                %%%%% 计算dissipative limiter Q
                Qx_i_plu_1_2=abs(ax_i_plu_1_2_RL)*(Cx_i_plu_1_2_R-Cx_i_plu_1_2_L);
                Qy_i_plu_1_2=abs(ay_i_plu_1_2_RL)*(Cy_i_plu_1_2_R-Cy_i_plu_1_2_L);

                Qx_i_min_1_2=abs(ax_i_min_1_2_RL)*(Cx_i_min_1_2_R-Cx_i_min_1_2_L);
                Qy_i_min_1_2=abs(ay_i_min_1_2_RL)*(Cy_i_min_1_2_R-Cy_i_min_1_2_L);
                
            case 'MTVDLF'
                Qx_i_plu_1_2=(max(abs(fc_x_i_plu_1_2_R),abs(fc_x_i_plu_1_2_L)))*(Cx_i_plu_1_2_R-Cx_i_plu_1_2_L);
                Qy_i_plu_1_2=(max(abs(fc_y_i_plu_1_2_R),abs(fc_y_i_plu_1_2_L)))*(Cy_i_plu_1_2_R-Cy_i_plu_1_2_L);
                Qx_i_min_1_2=(max(abs(fc_x_i_min_1_2_R),abs(fc_x_i_min_1_2_L)))*(Cx_i_min_1_2_R-Cx_i_min_1_2_L);
                Qy_i_min_1_2=(max(abs(fc_y_i_min_1_2_R),abs(fc_y_i_min_1_2_L)))*(Cy_i_min_1_2_R-Cy_i_min_1_2_L);
        end


        Fx_i_plu_1_2=0.5*(fc_x_i_plu_1_2_R+fc_x_i_plu_1_2_L-Qx_i_plu_1_2);
        Fy_i_plu_1_2=0.5*(fc_y_i_plu_1_2_R+fc_y_i_plu_1_2_L-Qy_i_plu_1_2);
        Fx_i_min_1_2=0.5*(fc_x_i_min_1_2_R+fc_x_i_min_1_2_L-Qx_i_min_1_2);
        Fy_i_min_1_2=0.5*(fc_y_i_min_1_2_R+fc_y_i_min_1_2_L-Qy_i_min_1_2);

        %%-------计算扩散�?------%%
        diffu_x=c1(i+1,j)-2*c1(i,j)+c1(i-1,j);
        diffu_y=c1(i,j+1)-2*c1(i,j)+c1(i,j-1);
        
        c2(i,j)=c1(i,j)-(dt/dx)*(Fx_i_plu_1_2-Fx_i_min_1_2)-(dt/dy)*(Fy_i_plu_1_2-Fy_i_min_1_2)+(dt/(dx^2))*D*diffu_x+(dt/(dy^2))*D*diffu_y;
       end
    end
    c1=c2;
end
toc;

%% 画图
contourf(x1,y1,c1')
% 设置坐标轴范�?
xlim([-3 3]);
ylim([-3 3]); 
disp(max(max(c1)));
disp(min(min(c1)));
saveas(gcf,  method + "+" + Limiter + '.svg', 'svg');
end

