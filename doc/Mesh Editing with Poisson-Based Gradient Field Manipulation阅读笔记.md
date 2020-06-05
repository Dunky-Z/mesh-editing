## Mesh Editing with Poisson-Based Gradient Field Manipulation阅读笔记
这篇论文介绍了一种基于泊松梯度场控制的网格编辑技术，运用该技术可以对三维网格进行变形，去噪，融合等操作。

说起泊松编辑，经典的应用应该是图像融合。图像融合的目的是为了将一张图片某一区域复制到目标图像中，并实现边缘的平滑过渡。

<center>
    <img src="https://gitee.com//dominic_z/markdown_picbed/raw/master/img/20200529104524.png"  width = "50% height = "50%">
</center>

### 图像泊松融合
<center>
    <img src=https://gitee.com//dominic_z/markdown_picbed/raw/master/img/20200529102214.png  width = "50% height = "50%">
</center>

图像泊松融合的主要思想是，根据原图像的梯度信息以及目标图像的边界信息，利用插值重新构建合成区域的像素。如上图，其中$u$表示原图需要合成进目标图像的部分，$V$是$u$的梯度场，$S$是合成后的目标图像，$\Omega$表示合成后被覆盖的区域，$\partial \Omega$表示其边界。设合并后的图像在$\Omega$区域的像素由未知函数$f$表示，以外的函数由已知函数$f^*$确定。  
图像合成的要求是合成后的图像看起来尽可能平滑的过渡，没有明显的边界，所以$\Omega$内的梯度应该尽可能的小。并且在边界上的像素值又和背景像素值相同，这就是一个带边界条件的二元泛函的极值问题，可以表述为下式：

$$\min _{f} \iint_{\Omega}\|\nabla f\|^{2}, \quad \text { s.t. }\left.f\right|_{\partial \Omega}=\left.f^{*}\right|_{\partial \Omega}$$

求解这样的极值问题需要使用变分法，本问题中的被积函数为：

$$F=\|\nabla f\|^2=f_x^2+f_y^2$$

带入欧拉-拉格朗日方程[<sup>[1]</sup>](#refer-anchor-1)（这是泛函取得极小值的必要条件的微分形式），则有：

$$\frac{\partial F}{\partial f}-\frac{\text{d}}{\text{d}x}(\frac{\partial F}{\partial f_x})-\frac{\text{d}}{\text{d}y}(\frac{\partial F}{\partial f_y})=0$$

因为$F$是$\nabla F$的函数，所以和$f$没关系，$\partial F/\partial f = 0$，此外：

$$\frac{\text{d}}{\text{d}x}(\frac{\partial{F}}{\partial{f_x}})=\frac{\text{d}}{\text{d}x}(2f_x)=2\frac{\partial^2f}{\partial{x^2}}$$

于是：

$$\frac{\partial^2{f}}{\partial x^2}+\frac{\partial^2{f}}{\partial{y^2}}=\Delta{f}=0$$

最后转化成了求解这个微分方程的。

以上的解法，虽然能够完成图像的融合，但是前景图片也会被进行一定的处理，导致前景图片与背景图片混合，不清晰。于是有了如下的处理方式：

$$\min _{f} \iint_{\Omega}\|\nabla f-V\|^{2}=\min _{f} \iint_{\Omega}\|\nabla f-\nabla u\|^{2}, \quad \text { s.t. }\left.f\right|_{\partial \Omega}=\left.f^{*}\right|_{\partial \Omega}$$

在原有条件下引入前景图片的梯度场作为引导场，为了使得合并后的图片中$\Omega$区域尽量接近前景图片$u$。合成后的图片在$\Omega$区域内的像素值$f$的梯度与$u$的梯度越接近，说明前景图片的原始纹理就保持得越好。此时被积函数为：

$$F=\|\nabla f-\nabla u\|^2=(f_x - u_x)^2+(f_y - u_y)^2$$

再次应用欧拉-拉格朗日方程，

$$\frac{\text{d}}{\text{d}x}[\frac{\partial{F}}{\partial{(f_x-u_x)^2}}]=\frac{\text{d}}{\text{d}x}[2(f_x-u_x)]=2(\frac{\partial^2f}{\partial{x^2}}-\frac{\partial^2{u}}{\partial{x^2}})$$

于是得到

$$\frac{\partial^{2} f}{\partial x^{2}}+\frac{\partial^{2} f}{\partial y^{2}}=\frac{\partial^{2} u}{\partial x^{2}}+\frac{\partial^{2} u}{\partial y^{2}}$$

对于一个向量场$F=P(x,y)\bm{i}+Q(x,y)\bm{j}$，它的散度记为$\operatorname{div} F$，且有：

$$\operatorname{div} F=\frac{\partial P}{\partial x}+\frac{\partial Q}{\partial y}$$

从散度的角度来定义拉普拉斯算子，拉普拉斯算子即是梯度的散度，即：

$$\Delta f=\operatorname{div}(\nabla f)=\nabla \cdot (\nabla f)=\nabla^2f$$

二维空间中，$f$是关于$x$和$y$的函数，所以：

$$\Delta f = \operatorname{div}(\frac{\partial{f}}{\partial x}\bm i,\frac{\partial{f}}{\partial{y}}\bm{j})=\frac{\partial^2{f}}{\partial{x^2}}+\frac{\partial^2{f}}{\partial{y^2}}$$

所以原式可以写成下列泊松方程的形式：

$$\Delta f=\operatorname{div}(\nabla u)$$


### 泊松网格编辑
回到三维的网格编辑问题，在二维空间，梯度和散度计算和表示都相对容易。如果我们能够找到在三维网格上梯度场，散度，拉普拉斯算子的定义，那泊松网格编辑的问题就迎刃而解了。

在三维网格编辑中仍然是求解这类似的泛函极值问题，

$$\min _{\phi} \iint_{\Omega}\|\nabla \phi-\mathbf{w}\|^{2} d A$$

只不过这里的引导场变成了定义在网格曲面上的分段常值矢量场。暂且将其也看做一个梯度场。这个引导矢量场可以通过顶点位置求得，所以可以当成是已知的。$\phi$是定义在网格上的一个函数。它与顶点息息相关。目的就是为了保证变形前后两个场的梯度尽可能保持一致就可以保证其细节特征不背破坏。只要通过改变顶点的位置，也就是$w$部分就可以通过前后梯度保持一致反推顶点的位置。


#### 三角网格的梯度
在每个三角形内部，其实可以通过三角形三个顶点的属性插值得到。如下所示，$f_i,f_j, f_k$分别是三个顶点关联的属性值，$f_u$可以通过$\alpha, \beta, \gamma$的线性插值得到，一般情况下，$\alpha, \beta, \gamma$分别可以表示为三个小三角形与大三角形的面积之比。


<center>
    <img src=https://gitee.com//dominic_z/markdown_picbed/raw/master/img/20191227144005445.png  width = "50% height = "50%">
</center>


因为$\alpha,\beta,\gamma$分别表示$v_i,v_j,v_k$的影响，因此将$\alpha,\beta,\gamma$定义为三个函数：$B_i(u),B_j(u),B_k(u)$,因为重心坐标的性质，$B_i(u)+B_j(u)+B_k(u)=1$.  
三角网格上任意一点都可以表示为：

$$f(u) = f_iB_i(u) + f_jB_j(u)+f_kB_k(u)$$

$f_i,f_j,f_k$分别为三个顶点的坐标信息。

梯度：

$$\nabla f(u) = f_i\nabla B_i(u) + f_j\nabla B_j(u)+f_k\nabla B_k(u)$$

因为$B_i(u)+B_j(u)+B_k(u)=1$，所以$\nabla B_i(u)+\nabla B_j(u)+\nabla B_k(u) = 0$

因此可得$\nabla B_i(u) = -\nabla B_j(u)-\nabla B_k(u)$

带入原式可得：

$$\nabla f(u)=(f_j-f_i)\nabla B_j(u)+(f_k-f_i)\nabla B_k(u)$$

问题来了如何求$B_j(u)$呢？下面以$\nabla B_i(u)$为例解释：

<center>
    <img src=https://gitee.com//dominic_z/markdown_picbed/raw/master/img/20191227120856240.png  width = "50% height = "50%">
</center>

如上左图，做一条平行于$v_j,v_k$的直线$l$，在平行线上取两个点：$v_x,v_y$，那么有：

$$f_x = f_iB_i(x)+f_jB_j(x)+f_kB_k(x)$$
$$f_y = f_iB_i(y)+f_jB_j(y)+f_kB_k(y)$$.

已知：

$$B_i(x) = s\triangle{xjk}/ s\triangle{ijk}$$ 

$$B_i(y)=s\triangle{yjk}/s\triangle{ijk}$$

因为$s\triangle{xjk}$,$s\triangle{yjk}$等底同高，所以$s\triangle{yjk} = s\triangle{xjk}$

因此：$B_i(x)=B_i(y)$

基于这个重要的结论，实际上所有平行于$v_x,v_y$的直线上的$B_i$均相等。为此，我们可以做出右侧的一个关于$B_i$的等值线。那么我们就可以知道，$B_i$梯度方向是垂直于等值线且指向顶点$v_i$的。在长度为底边高的线段上从0变化到1.那么，梯度的模长为高的倒数。因此：

$$\nabla B_{i}(u)=\left(x_{k}-x_{j}\right)^{\perp} /(2 A)$$

其中，$\perp$表示方向垂直底边指向另一个顶点。
在求梯度方向时可以顺便求出面法向，面法向的模的一半就是面积了，可以一起求出来。

将上式带入原式得：

$$\nabla f(u)=(f_j-f_i)(x_i-x_k)^{\perp}/(2A)+(f_k-f_i)(x_j-x_i)^{\perp}/(2A)$$


#### 三角面片顶点的散度
设向量值函数$w:S-R^3$,  $S$为网格模型。$w$表示各个顶点的向量，那么$w$的散度为，对顶点$v_i$ 1-邻域信息的求和  

$$\operatorname{div} \boldsymbol{w}\left(v_{i}\right)=\left.\sum_{T \in N\left(v_{i}\right)} \nabla B_{i}\right|_{T} \cdot \boldsymbol{w}_{T} \cdot A_{T}$$

其中$w_T$就是面片的梯度$\nabla f(u)$，$A_T$即面片面积。

因为三维空间，梯度有三个方向，所以散度需要将三个方向求和。


#### 顶点的拉普拉斯算子

$$\Delta f(v_i)=1/2\sum_{j\in{N_{v}(v_i)}}(\cot \alpha_{ij}+\cot \beta_{ij})(f_i-f_j)$$

#### 线性系统



$$\operatorname{Div}(\nabla \phi)=\operatorname{Div} \mathbf{w}$$

$$\mathbf A \mathbf f= \mathbf b$$


### 实验结果
#### 人体模型，固定脚部，头部上移

身高变高但是没有约束的部分有明显的变形

<center>
    <img src=https://gitee.com//dominic_z/markdown_picbed/raw/master/img/ori-left00.png width = "30% height = "30%" ><img src=https://gitee.com//dominic_z/markdown_picbed/raw/master/img/ori-left000.png width = "30% height = "30%">
</center>

#### Cubic模型,顶部拉伸
<center>
    <img src=https://gitee.com//dominic_z/markdown_picbed/raw/master/img/cubic-left010.png width = "30% height = "30%" ><img src=https://gitee.com//dominic_z/markdown_picbed/raw/master/img/cubic-left01.png width = "30% height = "30%">
</center>

#### Cylinder模型，顶部拉伸
<center>
    <img src=https://gitee.com//dominic_z/markdown_picbed/raw/master/img/Cylinder-down-00.png width = "30% height = "30%" ><img src=https://gitee.com//dominic_z/markdown_picbed/raw/master/img/Cylinder-down000.png width = "30% height = "30%">
</center>

<center>
    <img src=https://gitee.com//dominic_z/markdown_picbed/raw/master/img/Cylinder-left-01.png width = "30% height = "30%" ><img src=https://gitee.com//dominic_z/markdown_picbed/raw/master/img/Cylinder-left010.png width = "30% height = "30%">
</center>

 
### 附录
<div id="refer-anchor-1"></div>
数值方法-Chapter21.变分法初步