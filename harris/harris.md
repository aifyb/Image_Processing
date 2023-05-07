# Harris 角点检测算法

## 基本思想

图像上在某点沿着各个方向灰度值明显变化时，该点大概率被认为是角点。为了提高可靠性，从图像局部的小窗口观察图像特征，该窗口向任意方向的移动都导致图像灰度的明显变化，则认为此处为角点。

<img src="https://image.aifyb.tech/image_processing/image-20230402170955312.png" style="zoom:50%;" aligend="middle" />

## 数学表达

将图像窗口平移 $[u,v]$ 产生灰度变化 $E(u,v)$ ：
$$
E(u,v) =\sum\limits_{x,y}w(x,y)\left[ I(x+u,y+v)-I(x, y) \right]^2
$$
其中，$w(x,y)$ 为窗口函数，可以认为是各像素点灰度值对角点灰度值的贡献度，一般可取定义在窗口内的高斯函数或者常值函数。

二元函数 $f(x,y)$ 在点 $(x_0, y_0)$ 处的泰勒展开式如下给出，
$$
\begin{equation}
f(x,y)=f(x_0,y_0)+\frac{\partial f}{\partial x}(x_0,y_0)(x-x_0)+\frac{\partial f}{\partial y}(x_0,y_0)(y-y_0)+ \cdots
\end{equation}
$$
可以对 $I(x+u, y+v)$ 在点 $(x,y)$ 处进行泰勒展开，略去二次项及以上对其近似，则有
$$
\begin{equation}
\begin{split}
I(x+u,y+v) &= I(x,y)+u \frac{\partial f}{\partial x}(x,y)+v\frac{\partial f}{\partial y}(x,y)+ \cdots \\
&\approx I(x,y)+u \frac{\partial f}{\partial x}(x,y)+v\frac{\partial f}{\partial y}(x,y)\\
&=I(x,y)+u I_{x}(x,y)+v I_{y}(x,y)
\end{split}
\end{equation}
$$
$I_{x}(x,y), I_{y}(x,y)$ 可由图像与 Sobel 算子或 Prewitt 算子卷积得到，从而 $\left[ I(x+u,y+v)-I(x, y) \right]^2$ 可以化简为，
$$
\begin{equation}
\begin{split}
 \left[ I(x+u,y+v)-I(x, y) \right]^2 &\approx \left[ I(x,y)+u I_{x}(x,y)+v I_{y}(x,y)-I(x, y) \right]^2 \\
 &=\left[u I_{x}(x,y)+v I_{y}(x,y) \right]^2 \\
 &=u^2 I_{x}^2(x,y) +2uvI_{x}(x,y)I_{y}(x,y) + v^2 I_{y}^2(x,y)
\end{split}
\end{equation}
$$
省略自变量并改写为矩阵形式，
$$
\begin{equation}
    \begin{split}
     \left[ I(x+u,y+v)-I(x, y) \right]^2 &\approx u^2 I_{x}^2(x,y) +2uvI_{x}(x,y)I_{y}(x,y) + v^2 I_{y}^2(x,y) \\
     &=u^2 I_{x}^2 +2uvI_{x}I_{y} + v^2 I_{y}^2 \\
     &= \left( \begin{array}{cc} u, v  \end{array} \right) 
     \left( \begin{array}{cc}  
        I_x^2 &   I_xI_y \\
        I_xI_y & I_y^2  
    \end{array} \right) \left(\begin{array}{c} u \\v \end{array}\right)
    \end{split} 
\end{equation}
$$
将上式代入灰度变化函数 $E(u,v)$ 有，
$$
\begin{equation}
    \begin{split}
    E(u,v) &=\sum\limits_{x,y}w(x,y)\left[ I(x+u,y+v)-I(x, y) \right]^2 \\
    &\approx\sum\limits_{x,y}w(x,y) \left( \begin{array}{cc} u, v  \end{array} \right) 
     \left( \begin{array}{cc}  
        I_x^2 &   I_xI_y \\
        I_xI_y & I_y^2  
    \end{array} \right) \left(\begin{array}{c} u \\v \end{array}\right) \\
    &=\left( \begin{array}{cc} u, v  \end{array} \right)\sum\limits_{x,y}w(x,y)  \left( \begin{array}{cc}  
        I_x^2 &   I_xI_y \\
        I_xI_y & I_y^2  
    \end{array} \right) \left(\begin{array}{c} u \\v \end{array}\right) \\
    &=\left( \begin{array}{cc} u, v  \end{array} \right)M \left(\begin{array}{c} u \\v \end{array}\right)
    \end{split} 
\end{equation}
$$
其中，$M=\sum\limits_{x,y}w(x,y)  \left( \begin{array}{cc}  
        I_x^2 &   I_xI_y \\
        I_xI_y & I_y^2  
    \end{array} \right)$ ，称为海森矩阵，显然 $M$ 为实对称矩阵，根据实对称矩阵必可相似对角化的性质我们有，
$$
M=P\Lambda P^{-1} = P\Lambda P^{T}=P \left( \begin{array}{cc}  
        \lambda_1 &   0 \\
        0 & \lambda_2  
    \end{array} \right)P^{T}
$$
其中，$\lambda_1,\lambda_2$ 为 $M$ 的两个特征值，$P$ 为正交矩阵，我们知道，正交矩阵的物理意义类似于旋转矩阵，即正交矩阵作用在一个向量上相当于将该向量做了一个旋转（这部分内容有时间我再单开一篇文章介绍）。

进而将 $E(u,v)$ 继续改写为，
$$
\begin{equation}
\begin{split}
    E(u,v) &\approx \left( \begin{array}{cc} u, v  \end{array} \right)M \left(\begin{array}{c} u \\v \end{array}\right)\\
&=\left( \begin{array}{cc} u, v  \end{array} \right)P \left( \begin{array}{cc}  
        \lambda_1 &   0 \\
        0 & \lambda_2  
    \end{array} \right)P^{T}\left(\begin{array}{c} u \\v \end{array}\right)\\
&= \left( \begin{array}{cc} u, v  \end{array} \right)P \left( \begin{array}{cc}  
        \lambda_1 &   0 \\
        0 & \lambda_2  
    \end{array} \right) \left(\left( \begin{array}{cc} u, v  \end{array} \right)P\right)^T\\
&=\left( \begin{array}{cc} u', v'  \end{array} \right) \left( \begin{array}{cc}  
        \lambda_1 &   0 \\
        0 & \lambda_2  
    \end{array} \right) \left(\begin{array}{c} u' \\v' \end{array}\right)\\
    &=\lambda_1u'^2 + \lambda_2 v'^2 \\
    &=\frac{u'^2}{\frac{1}{\lambda_1}} + \frac{v'^2}{\frac{1}{\lambda_2}}
\end{split} 
\end{equation}
$$
很明显，当灰度值变化量为给定的常数c即 $E(u,v)=c$ 时，$\frac{u'^2}{\frac{1}{\lambda_1}} + \frac{v'^2}{\frac{1}{\lambda_2}}=c$ 为以 $(u', v')$ 为变量的椭圆，从而说明 $E(u,v)$ 为椭圆区域，可以看作将 $\frac{u'^2}{\frac{1}{\lambda_1}} + \frac{v'^2}{\frac{1}{\lambda_2}}$ 旋转得到，其半轴的长度分别为 $\sqrt{\frac{1}{\lambda_1}}=\lambda_1^{-\frac{1}{2}}$ 和 $\sqrt{\frac{1}{\lambda_2}}=\lambda_2^{-\frac{1}{2}}$ 。短轴地方，代表像素值变化最剧烈的方向，因为得到相同灰度所需偏移量小；长轴地代表像素值变化最缓慢的地方，因为得到相同灰度所需偏移量大。

椭圆的长短轴分3种情况：

- $\lambda_1$ 和 $\lambda_2$ 都很小，即椭圆的长短轴都很大且相近，代表所有方向的像素值变化都很缓慢，说明窗口在平坦区域；

- $\lambda_1$ 和 $\lambda_2$ 一个大一个小，即椭圆的长轴很长，短轴很短，代表在单一方向像素值变化快，其他方向变化缓慢，说明窗口包含边缘；

- $\lambda_1$ 和 $\lambda_2$ 都很大，即椭圆的长短轴都很小且相近，代表两个（多个）方向的像素值变化都很缓快，说明窗口内包含角点；

  <img src="https://image.aifyb.tech/image_processing/image-20230402200423893.png" aligen="middle"/>

为了体现 $\lambda_1$ 和 $\lambda_2$ 都很大，从而判断窗口内是否有角点，我们定义角点响应函数 $R$ ：
$$
\begin{equation}
R \xlongequal{def} \det (M) -k( \operatorname{trace}(M) )^2
\end{equation}
$$
 其中，$\det(M) = \lambda_1 \lambda_2$ , $\operatorname{trace}(M)=\lambda_1 +\lambda_2$ , k为经验常数，一般取 $0.04\sim 0.06$ 。

- 增大k值，降低角点检测的灵敏度，减少被检角点的数量；
- 减少k值，增加角点检测的灵敏度，增加被检测角点的数量；
- R只与M的特征值有关
  - 角点：R为大数值正数
  - 边缘：R为大数值负数
  - 平坦区域：R为小数值（可正可负）

## Harris 角点检测算法

具体流程：

1. 彩色图像转换为灰度图像，加快处理速度；
2. 利用 sobel 算子计算整张图的 $I_x, I_y$ ；
3. 构建海森矩阵，将高斯滤波器分别作用于 $I^2_x, I_xI_y, I^2_y$ ；
4. 计算响应，计算每个像素的R；
5. 非极大值抑制。

## 总结与补充

### Harris角点的性质

旋转不变性：椭圆旋转过一定角度但是其形状保持不变，即特征值不变；

- 焦点响应函数R对于图像的旋转具有不变性；
- 对于图像灰度的仿射变化具有部分不变性；$I\rightarrow aI+b$ 
- 对于图像几何尺度变化具有不变性；随尺度变化，Harris 角点检测的性能下降；

### SHI-TOMASI 角点

由于 Harris角点检测算法的稳定性和k值有关，而 k是个经验值，不好设定最佳值，Shi-Tomasi发现，角点的稳定性和矩阵 M的较小特征值有关，于是直接用较小的那个特征值作为分数，这样就不需调整 k值。
$$
R = \min(\lambda_1, \lambda_2)
$$



