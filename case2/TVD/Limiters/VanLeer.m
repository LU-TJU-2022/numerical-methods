function [Limiter1x, Limiter1y] = VanLeer(c1, dx, dy, i, j)
    % VanLeer 限制器计算函数
    % 输入:
    %   c1  - 计算网格上的浓度场矩阵
    %   dx, dy - 网格间距
    %   i, j - 当前计算点的索引
    % 输出:
    %   Limiter1x - x 方向的 VanLeer 限制器值
    %   Limiter1y - y 方向的 VanLeer 限制器值
    
    XHx_i=(c1(i,j)-c1(i-1,j))/(c1(i+1,j)-c1(i,j));%%%seita(i)
    XHy_i=(c1(i,j)-c1(i,j-1))/(c1(i,j+1)-c1(i,j));%%%seita(i)
    %%%%%%%计算limiter——i-1
    XHx_i_min_1=(c1(i-1,j)-c1(i-2,j))/(c1(i,j)-c1(i-1,j));%%%seita(i-1)
    XHy_i_min_1=(c1(i,j-1)-c1(i,j-2))/(c1(i,j)-c1(i,j-1));%%%seita(i-1)
    %%%%%%%%%计算limiter——i+1
    XHx_i_plu_1=(c1(i+1,j)-c1(i,j))/(c1(i+2,j)-c1(i+1,j));%%%seita(i+1)
    XHy_i_plu_1=(c1(i,j+1)-c1(i,j))/(c1(i,j+2)-c1(i,j+1));%%%seita(i+1)
    
     if isnan(XHx_i)  %%%  XHx = 0/0
        Limiter1x_i=2*(c1(i+1,j)-c1(i,j))/dx;
     elseif XHx_i<0     %%%  XHx<0  or  %%%XHx = -inf
        Limiter1x_i=0;
     else    %%%  XHx>0  or   XHx_i = inf
        Limiter1x_i=(2-2/(1+XHx_i))*(c1(i+1,j)-c1(i,j))/dx;
     end
     if isnan(XHy_i)
         Limiter1y_i=2*(c1(i,j+1)-c1(i,j))/dy;
     elseif XHy_i<0
        Limiter1y_i=0;
     else
        Limiter1y_i=(2-2/(1+XHy_i))*(c1(i,j+1)-c1(i,j))/dy;
     end

     if isnan(XHx_i_min_1)
        Limiter1x_i_min_1=2*(c1(i,j)-c1(i-1,j))/dx;%%%sigama(i-1)
     elseif XHx_i_min_1<0
        Limiter1x_i_min_1=0;
     else
        Limiter1x_i_min_1=(2-2/(1+XHx_i_min_1))*(c1(i,j)-c1(i-1,j))/dx;%%%sigama(i-1)
     end
     if isnan(XHy_i_min_1)
        Limiter1y_i_min_1=2*(c1(i,j)-c1(i,j-1))/dy;%%%sigama(i-1)
     elseif XHy_i_min_1<0
        Limiter1y_i_min_1=0;
     else
        Limiter1y_i_min_1=(2-2/(1+XHy_i_min_1))*(c1(i,j)-c1(i,j-1))/dy;%%%sigama(i-1)
     end

     if isnan(XHx_i_plu_1)
        Limiter1x_i_plu_1=2*(c1(i+2,j)-c1(i+1,j))/dx;%%%sigama(i+1)
     elseif XHx_i_plu_1<0
        Limiter1x_i_plu_1=0;
     else
        Limiter1x_i_plu_1=(2-2/(1+XHx_i_plu_1))*(c1(i+2,j)-c1(i+1,j))/dx;%%%sigama(i+1)
     end
     if isnan(XHy_i_plu_1)
        Limiter1y_i_plu_1=2*(c1(i,j+2)-c1(i,j+1))/dy;%%%sigama(i+1)
     elseif XHy_i_plu_1<0
        Limiter1y_i_plu_1=0;
     else
        Limiter1y_i_plu_1=(2-2/(1+XHy_i_plu_1))*(c1(i,j+2)-c1(i,j+1))/dy;%%%sigama(i+1)
     end

    % 返回结果
    Limiter1x = [Limiter1x_i; Limiter1x_i_min_1; Limiter1x_i_plu_1];
    Limiter1y = [Limiter1y_i; Limiter1y_i_min_1; Limiter1y_i_plu_1];
end
