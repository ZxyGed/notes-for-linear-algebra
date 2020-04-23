# MIT线性代数

课程：[麻省理工线性代数](https://www.bilibili.com/video/BV1Kb411M72F)

## 方程组的几何解释

从行的角度看，矩阵可以看作是方程组（方程组也可以看成面，于是化为面的交点）；

从列的角度，每一列可以看作一个三维向量，化为三个向量的线性和，于是对于如下问题

:star: `Can I solve Ax = b for every b ? ` $（A_{m\times n}）$ $\Longleftrightarrow$ `Do the linear combinations of the columns fill the M-D space?` 这需要从列空间理解。

换句话说，想要 `Ax = b` 有解，b需要属于A的列空间，即 $b \in C(A)$

当三个向量共面时，便无法充满三维空间，这种情况成为奇异（singular case），此时矩阵不可逆

:star:矩阵乘法的核心是向量的线性组合！矩阵在左进行行变换，矩阵在右进行列变换

-  $matrix\times column=column$，此时从列的角度可以看作列向量的线性组合

$$
\begin{bmatrix} column1 & column2 & column3\end{bmatrix}
\cdot
\begin{bmatrix} x\\ y\\z\end{bmatrix}=
x \cdot column1 + y \cdot column2 + z \cdot column3
$$

- $row \times matrix = row$，此时从行的角度也可以看作行向量的线性组合

$$
\begin{bmatrix} x & y & z\end{bmatrix}
\cdot
\begin{bmatrix} row1\\row2\\row3\end{bmatrix}
=
x \cdot row1 + y \cdot row2 + z \cdot row3
$$



## 矩阵消元

利用矩阵乘法是线性组合的特性，对于 $row\times matrix$ 的情况，可以将其进行行角度的拆分
$$
\begin{bmatrix} operation1 \\ operation2 \\ operation3 \end{bmatrix} 
\cdot 
\begin{bmatrix} row1 \\ row2 \\ row3 \end{bmatrix} = 
\begin{bmatrix} linear\_combins1 \\ linear\_combins2 \\ linear\_combins3 \end{bmatrix}
$$
即，基于 $row \times matrix = row$ ，将matrix变为三次线性组合结果的叠加，对于第一行，如果 $x1=1,y1=0, z1=0$ 则线性组合输出第一行结果为原第一行，如果 $x2=0,y2=1, z2=0$ 则第二行保持不变，第三行同理，由此便退出了单位矩阵的概念，即对自身操作后仍未自身。

**标记语言**：E~21~表示第二行与第一行之间的线性变化（E表示初等矩阵element），同理表示E~32~，于是整体的计算过程可以表示为 $E_{32}(E_{21}A)=U$ 其中A为系数矩阵，U为上三角矩阵（upper）

通过结合律，我们可以得到更精简的形式 $E_{32}(E_{21}A)=(E_{32}E_{21})A=U$



## 乘法和逆矩阵

对于 $A\cdot B=C$ ，具有如下四种理解：

1. 从元素角度进行理解，其满足

   $c_{ij}=row_i\times column_j=a_{i1}b_{1j}+\dots+a_{in}b_{nj}=\sum_{k=1}^na_{ik}b_{kj}$

2. 从列的角度，可以看作右边矩阵的每一列对左边进行线性组合

3. 从行的角度，可以看作左边矩阵的每一行对右边进行线性组合

4. 之前都是 $row \times column$ 的形式，转换思路为 $column \times row$ ，这种形式仍可以从行列的角度去理解，这就转换为行空间与列空间，即从列的角度，运算结果都是列向量的倍数，从行的角度，运算结果都是行向量的倍数，均在一条直线上



**矩阵运算技巧：分块乘法**
$$
\begin{bmatrix} A1&A2\\A3&A4 \end{bmatrix} \cdot
\begin{bmatrix} B1&B2\\B3&B4 \end{bmatrix}=矩阵运算结果，因为每一块仍可以看作一个线性组合
$$


**逆矩阵的理解**

逆矩阵的产生是对于 $EA=U$ 想要寻找一个方法使其 $E^{-1}EA=U \Leftrightarrow E^{-1}E=I$

观察如下矩阵,会发现左边两个矩阵互为逆矩阵，先理解中间矩阵的含义，即一三行不变，第二行减去三倍的第一行，如果要将第二行变回原样，只需重新加上三倍的第一行即可，这就是其逆矩阵做的事，由此便可以很好地理解逆矩阵的含义
$$
\begin{bmatrix}1&0&0 \\ 3&1&0 \\ 0&0&1\end{bmatrix}
\cdot
\begin{bmatrix} 1&0&0 \\ -3&1&0 \\ 0&0&1\end{bmatrix}
=\begin{bmatrix}1&0&0 \\ 0&1&0 \\ 0&0&1\end{bmatrix}
$$



**奇异不可逆的矩阵（singular case, no inverse）**：`you can find a vector with Ax=0` 即$Ax=0\Leftrightarrow A^{-1}Ax=0\rightarrow x=0$ 矛盾，因此对于该类型的矩阵不可逆



**Gauss-Jordan elimination**

$E\cdot \begin{bmatrix} A& I\end{bmatrix}=\begin{bmatrix} I&?\end{bmatrix}\qquad\because E\cdot A=I\quad \therefore E=A^{-1}\ \rightarrow\ \ ?\ is\ A^{-1}$ 

$E\cdot \begin{bmatrix} A& I\end{bmatrix}=\begin{bmatrix} I&A^{-1}\end{bmatrix}$ 即将增广矩阵添加到右侧将系数矩阵转换为单位矩阵即得逆矩阵



## A的LU分解

即 $A=LU\rightarrow L=E^{-1}$ ,L代表lower，下三角矩阵，可以这样理解，在将A化成U的过程中，第n行减去第n-1行，也就是最开始的操作是不会涉及n+1之后的行的，因此E本身就是下三角矩阵，其逆同样也是下三角矩阵。

$AA^{-1}=I$ 将两边同时转置（左边顺序会发生颠倒），得到$(A^{-1})^TA^T=I\rightarrow (A^{-1})^T=(A^T)^{-1}$ 即对于单个矩阵而言，其转置和逆运算可以颠倒

LU分解相对于高斯消元的优势在于对于大矩阵，在b变化的时候也无需全部重新计算，因为计算过程分为求逆以及回代，复杂度分别为 $n^3$ 以及 $n^2$ ，高斯方法在b改变情况下都需要重算，LU方法仅需要重新回代，具体意义参考 [线性代数笔记10——矩阵的LU分解](https://blog.csdn.net/sunbobosun56801/article/details/82190610)

## 转置-置换-向量空间

公式：$(AB)^T=B^TA^T$

对于n维矩阵，共有 n! 个置换矩阵P，其中 $P^{-1}=P$，同时由于$P^T P=I\rightarrow P^T=P^{-1}$

对称矩阵：$P=P^T$ 事实上对于矩阵R，$RR^T$ 结果永远都是对称矩阵（可以通过计算过程理解），可以将其表示为 $RR^T=(RR^T)^T=(R^T)^TR^T=RR^T$

向量空间：满足其中的向量相乘相加后仍旧处于其中（对加法和乘法是封闭的，即线性组合封闭），这也就意味着所有向量空间毕竟包含零向量，如$R^n$ 表示由n维向量构成的向量空间 

向量子空间（列空间C(A),C 代表column）：也就是向量空间的一部分，但其满足线性组合封闭的性质，例如对于一位的向量空间而言，其向量子空间有零向量空间，本身（过原点 的直线）；对于二位向量空间，其子空间包括零向量空间，过原点的直线以及本身；对于三维向量空间，其子空间包括零向量空间，过原点的直线，过原点的平面以及自身。。。以此类推

:star:零向量空间指的是 $\{x|Ax=0 \}$，其必定是向量空间，首先其包含零向量，其次，对于 $w,v \in x$，其满足线性组合封闭性，而对于 $\{x|Ax=b,b\neq0\}$ 必定不是向量空间，因为其不包含零向量，对于零向量空间的表示，通过取特殊解可以理解（可能为过原点的平面或过原点的直线）

考虑矩阵 $A=\begin{bmatrix}1&3\\2&3\\4&1\end{bmatrix}$ ，其列（两个三维向量）构成列空间（三维向量空间中过原点的平面），同理对于 $A_{10\times 5}$ 其构成的是十维向量空间中，过原点的五维列空间（需要确保每列线性无关）

## 列空间和零向量

以 $R^3$ 为例，考虑其中的子空间P(plane)以及L(line)，$P\cup L$ 并不构成向量子空间，因为他们不满足加法封闭，而$P\cap L$ 是满足的，考虑 $w,v \in P\cap L$，则 $w+v\in\ P且\ w+v\in\ L$  因为他们同时属于子空间P与L中，乘法同理，因此其交集属于向量子空间

其余内容结合上一节理解

## 求解 $Ax=0$ 主变量、特解

rref(reduce row echelon form，即简化行阶梯矩阵形式)：[Python3 矩阵求最简行阶梯矩阵](https://blog.csdn.net/ikoiiii/article/details/92075982?depth_1-utm_source=distribute.pc_relevant.none-task-blog-BlogCommendFromBaidu-3&utm_source=distribute.pc_relevant.none-task-blog-BlogCommendFromBaidu-3)

其中U为行阶梯矩阵，R为最简行阶梯矩阵
$$
A=\begin{bmatrix}
1&2&2&2\\
2&4&6&8\\
3&6&8&10
\end{bmatrix}
\rightarrow U=\begin{bmatrix}
1&2&2&2\\
0&0&2&4\\
0&0&0&0
\end{bmatrix}
\rightarrow R=\begin{bmatrix}
1&2&0&-2\\
0&0&1&2\\
0&0&0&0
\end{bmatrix}
$$
矩阵的**秩**为每个矩阵的非零行数目，其中非零行的首个非零元素对应的变量称为**主元**，其余变量成为**自由变量**，在上式中$x_1$和$x_3$是主元，$x_2$和$x_4$是自由变量

可以依次取自由变量中单独一个量为1，其余为0来计算特解，而最终的零向量空间为其所有特解的线性组合，如上式中$x_2$和$x_4$是自由变量，则令解分别为$\begin{bmatrix}?\\0\\?\\1 \end{bmatrix}$ 以及 $\begin{bmatrix}?\\1\\?\\0 \end{bmatrix}$，求出特解，最终解为$x=c\begin{bmatrix}-2\\1\\0\\0 \end{bmatrix}+d\begin{bmatrix}2\\0\\-2\\1 \end{bmatrix}$ 

:star2: insight：将R中主列以及自由列重新排列，可以得到如下形式
$$
R=\begin{bmatrix}
1&0&2&-2\\
0&1&0&2\\
0&0&0&0
\end{bmatrix}
=\begin{bmatrix}
I&F\\
0&0
\end{bmatrix}（F\ means\ free\ column）
$$

$$
Rx=0\rightarrow\begin{bmatrix}I&F\\0&0 \end{bmatrix}x=0\quad \Rightarrow \quad x=\begin{bmatrix}-F\\ I \end{bmatrix}=\begin{bmatrix}-2&2\\0&-2\\1&0\\0&1\end{bmatrix}
$$

$x$ 的列空间（零向量空间），为特解组成的矩阵，对比原答案，发现满足，事实上，为了求解方便，无需将R进行列重排，在原基础上$\begin{bmatrix}?&?\\0&1\\?&?\\1&0 \end{bmatrix}$ 将 -F 依次填入即可

