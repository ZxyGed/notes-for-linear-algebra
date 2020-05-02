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



**奇异不可逆的矩阵（singular case, no inverse）**：`you can find a vector with Ax=0` 即$Ax=0\Leftrightarrow A^{-1}Ax=0\rightarrow x=0$ 矛盾，因此对于该类型的矩阵不可逆，即 $A$ 零空间只有零向量，即 $A$ 各列线性无关



**Gauss-Jordan elimination**

$E\cdot \begin{bmatrix} A& I\end{bmatrix}=\begin{bmatrix} I&?\end{bmatrix}\qquad\because E\cdot A=I\quad \therefore E=A^{-1}\ \rightarrow\ \ ?\ is\ A^{-1}$ 

$E\cdot \begin{bmatrix} A& I\end{bmatrix}=\begin{bmatrix} I&A^{-1}\end{bmatrix}$ 即将增广矩阵添加到右侧将系数矩阵转换为单位矩阵即得逆矩阵



## A的LU分解

即 $A=LU\rightarrow L=E^{-1}$ ,L代表lower，下三角矩阵，可以这样理解，在将A化成U的过程中，第n行减去第n-1行，也就是最开始的操作是不会涉及n+1之后的行的，因此E本身就是下三角矩阵，其逆同样也是下三角矩阵。

$(AB)^T=B^TA^T$ 可以从行列变换角度理解，即矩阵在左进行的是行变换，当进行转置，行变为列，于是需要位置互换，即行变换转变为列变换。

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

## 求解Ax=b可解性以及解的结构

$Ax=b\Rightarrow Ax=b+0$ 即其解为一个自由变量都为0的特解以及零向量空间组成，从空间上看，是将0向量空间平移到正好穿过特解，特解为particular（前提是特解存在）

:star2:：以下讨论秩r与解的情况

首先需要明白对于 $A_{m\times n}$，m为方程式个数，n为变量个数，$r\leq min(m,n)$

当 $r<n$，存在自由变量，零向量空间包含无穷多向量，否则仅包含全为0的向量

当 $r<m$，存在全为0的行，这意味着存在无特殊解的情况

| 条件      | 最简阶梯矩阵                             | 解的情况         |
| --------- | ---------------------------------------- | ---------------- |
| $r=m=n$   | $R=I$                                    | 有且仅有一解     |
| $r=n<m$   | $R=\begin{bmatrix}I\\0 \end{bmatrix}$    | 无解或者仅有一解 |
| $r=m<n$   | $R=\begin{bmatrix} I&F\end{bmatrix}$     | 有无穷解         |
| $r<m,r<n$ | $R=\begin{bmatrix}I&F\\0&0\end{bmatrix}$ | 无解或无穷解     |

1. $r=m=n$，此时矩阵满秩，$r=n$ 推出该矩阵覆盖n维空间同时其不存在自由变量，零向量空间仅包含全为零的向量，又由于 $m=n$，则 $R=I$ ，不存在全为0的行，即一定有1解，综合来看，在该情况下有且仅有一解
2. $r=n<m$，此时无自由变量，零向量空间仅包含全为零的向量，又$R=\begin{bmatrix}I\\0 \end{bmatrix}$，因此取决于特殊解是否存在，综上无解或者仅有一解。
3. $r=m<n$，此时存在自由变量，即一定有无穷解（零向量空间），$R=\begin{bmatrix} I&F\end{bmatrix}$，一定存在特殊解，综上有无穷解
4. $r<m,r<n$，此时存在自由变量，同时存在全为0的行，$R=\begin{bmatrix}I&F\\0&0\end{bmatrix}$ ，无解解或无穷解

##  线性相关性、基、维数

**线性无关**：如果找不到一个线性组合，使其对向量 $x_1,x_2,\dots,x_n$ 而言，使得 $c_1x_1+c_2x_2+\dots+c_nx_n=0$，则这些向量线性无关的（除去所有系数为0的情况）。

:star:换句话说，求解 $Ax=0,A=[x_1,x_2,\dots,x_n]$ ，如果其零向量空间中包含非零向量，则其线性相关 $\Longleftrightarrow$ **如果 $r<n$ ，其必定线性相关**

重述：when $v_1,v_2,\dots,v_n$ are columns of A, They are independent if null space of A is {zero vector} They are dependent if $Ac=0$ for some none-zero vector c in the null space.

:star:由定义可知，只要包含零向量，则其必定线性相关，因为零向量与任意向量线性相关，即取0向量的系数非零，其他向量系数为0

**生成空间定义**：Vectors $v_1,v_2,\dots ,v_l$ span a space means the space consists of all combinations of those vectors.

**基(basis)定义**：Basis for a space is a sequence of vectors $v_1,v_2,\dots,v_n$ with 2 properties:

1.  they are independent
2. they span the space

考虑 $\begin{bmatrix}1\\1\\2 \end{bmatrix}$ 以及 $\begin{bmatrix}2\\3\\4 \end{bmatrix}$ 也能成为一组基，其为过这两个向量的平面的基

:star:对于 $R^n$ 中的n个向量，想要组成一组基，由这n个向量为列向量组成的矩阵必须可逆（因为可逆的判断条件同线性无关的判断条件相同，参见乘法和逆矩阵那一节）

**维度的定义**：Given a space,every basis for the space has the number of vectors, the number is dimension of the space.

:star:rank(A) = the number of pivot columns = dimension of C(A) [列空间]

对应的，N(A)[零向量空间]的维度为自由变量的个数

## 四个基本子空间

[10. MIT线性代数---四个基本子空间](https://zhuanlan.zhihu.com/p/44313005)

:star:行变换的过程不影响方程的解，行空间不变（因为做的是线性变换），但是列空间会变。

对于矩阵 $A_{m\times n}$

| 所处向量空间 | 空间                                | 维度  |
| :----------: | :---------------------------------- | :---: |
|    $R^m$     | column space：$C(A)$                |  $r$  |
|    $R^m$     | null space of A transpose：$N(A^T)$ | $m-r$ |
|    $R^n$     | row space：$C(A^T)$                 |  $r$  |
|    $R^n$     | null space：$N(A)$                  | $n-r$ |

注意左零空间（$N(A^T)$）的基的求解方法：为了方便求解（而不是转置重新计算），根据零空间的定义如下转换，$A^Ty=0\rightarrow y^TA=0$，此时y由行向量变为列向量。

此时使用Gauss-Jordan消元法，复刻出 $A\rightarrow R$ 的变化矩阵
$$
\big[A_{m\times n}\ |\ I_{m\times n}\big]\Longrightarrow \big[R_{m\times n}\ |\ E_{m\times n}\big]
$$
$E_{m\times n}$ 即为所求，满足 $EA=R$，以如下具体例子为例
$$
EA=\begin{bmatrix}-1&2&0\\1&-1&0\\-1&0&1\end{bmatrix}\cdot 
\begin{bmatrix}1&2&3&1\\1&1&2&1\\1&2&3&1\end{bmatrix}
=\begin{bmatrix}1&0&1&1\\0&1&1&0\\0&0&0&0\end{bmatrix}
=R
$$
首先可以知道 $r=2$，$N(A^T)$ 的维度为 $m-r=1$， 观察最后一行 $[-1\quad0\quad 1]$ ，其满足 $y^TA=0$，其即为左零空间的基

**矩阵空间**：将 $R^n$ 扩展为 $R^{n\times n}$ ，因为矩阵同时也满足向量空间的运算规则，所以可以把矩阵也看做空间，其中包括上三角矩阵、对称矩阵、对角矩阵等（大小一次减小）

## 矩阵空间、秩1矩阵和小世界图

[11. MIT线性代数---矩阵空间、秩1矩阵和小世界图](https://zhuanlan.zhihu.com/p/44500497)

$3\times 3$ 矩阵的维度为9，即在每个位置依次取1

上三角矩阵 (U) 以及对称矩阵 (S) 的维度为6，对角矩阵矩阵 ($D=S\cap U$) 的维度为3

得出公式 $dim(S)+dim(U)=dim(S\cap U)+dim(S+U)$ 即 $6+6=3+9$ ，注意此处是 $S+U$ 而非 $S \cup U$

微分方程与解空间：以 $\frac{d^2y}{dx^2}+y=0$ 为例，$y=c_1cosx+c_2sinx$，$cosx和sinx$ 是其一组基，解空间的维度是2

关于秩1矩阵、左零空间以及小世界图的内容，参见上述链接，很详细

## 图和网络

[12. MIT线性代数---图和网络](https://zhuanlan.zhihu.com/p/45189924)

从电势的理解，$x$ 表示每点的电势，$A$ 表示关联矩阵（流出为-1，流入为1）

$Ax$ 表示每条边的电势差，而 $A^T$ 表示每个节点的流入流出情况，在平衡情况下满足基尔霍夫定律，即流入流出平衡，即 $A^Tx=0$

从图的角度而言，存在环路的三条边线性相关，而完全线性无关的边得出的是一棵树（因此秩为节点数-1），此时每添加一条边，便添加一个环路 $dim(N(A^T))\sim \#free\ edges\sim \#loops$ ，左零空间的维数就是相互无关的回路的数量

考虑 $A^Tx=0$ 的左零空间，即 $dim(N(A^T))=m-r$ 

转化为图论语言可以理解为 $\#loops=\#edges-(\#nodes-1)$ ，移项得 $\#nodes-\#edges+\#loops=1$ ，即欧拉公式，从维度角度看，$零维-一维+二维=1$

整合一下三个公式：

- $e=Ax$ ，$e$为电势差矩阵，$A$为关联矩阵，$x$为每点的电势
- $f=Ce$ ，来源欧姆定律，即电势差与电流的关系，$e$为电势差，$f$为电流
- $A^Tx=0$ ，基尔霍夫公式，$x$表示电流，$A^T$ 为关联矩阵的转置

考虑外加电源的情况下，可以理解为外部电流的流入，转化为 $A^Tx=b$

得到最终平衡状态下的公式 $A^TCAx=b$

## 复习一

[13. MIT线性代数---复习一](https://zhuanlan.zhihu.com/p/45192869)

强烈建议所有题目都过一遍

交换矩阵的某两行，零空间不变是因为，零空间是解空间，交换某两行，对解空间没有影响。

## 正交向量与子空间

[14. MIT线性代数---正交向量与向量子空间](https://zhuanlan.zhihu.com/p/45193986)

正交：$x^Ty=0$ 

勾股定理：$\parallel \vec x\parallel^2+\parallel \vec y\parallel^2=\parallel \vec x+\vec y\parallel^2$

证明：$\parallel x+y\parallel^2=(\vec x+\vec y)^T(\vec x+\vec y)=x^Tx+y^Ty+x^Ty+xy^T=\parallel \vec x\parallel^2+\parallel \vec y\parallel^2$

**子空间的正交**：两个子空间正交，两个空间中的所有向量正交，因此对于两个相交的子空间，其必定不是正交的，因为对于空间的交线，已经不满足对于所有空间中向量正交。

**行空间与零空间正交，二者互相为n维空间的正交补（因为$rank(C(A^T))+rank(N(A))=n$)，列空间与左零空间正交，二者互为n空间的正交补**

因此在三维空间中，一组正交直线无法成为行空间及零空间，因为 $1+1\neq3$ ，二者不互补

**求解无解方程 $Ax=b$** 

为了清除坏数据，常常将方程组化为 $A^TA\cdot \hat x=A^Tb$ ，即投影到最近的位置（$A^T(b-A\hat x)=0$，详见下一节），给出如下结论：

$N(A^TA)=N(A),rank(A^TA)=rank(A)$ 

## 子空间的投影

[15. MIT线性代数---投影](https://zhuanlan.zhihu.com/p/45246414)

这篇文章已经总结的很好了，主要补充两点

1. P是投影矩阵，对于 $P^2=P$ 可以直接从图形上理解，即一次投影之后已经在列空间中，再次投影对结果没有影响
2. 对于不能进行如下化简 $P^2=A(A^TA)^{-1}A^T=AA^{-1}(A^T)^{-1}A^T=I$ ，主要原因在于A不是方阵，此时 $A^{-1}$ 不存在，在 $A$ 为可逆方阵时，列空间充满n维空间，此时对于任意b，均在列空间中



重要公式推导：

对于在线上的投影

$a^T(b-p)=a^T(b-ax)=0\rightarrow x=\frac{a^Tb}{a^Ta}$ 则 $p=ax=\frac{aa^T}{a^Ta}b=Pb$ 表示 $b$ 在 $a$ 上的投影，$P$ 为投影矩阵，其中通过维数计算可以判断 x 为 $n\times 1$ 的向量

对于在平面上的投影

$A^T(b-Ax)=0\rightarrow A^Tb=AA^Tx\rightarrow x=(A^TA)^{-1}A^Tb$ 此处不直接除是由于 $A^TA$ 是一个可逆方阵，之前在线上的投影直接除是因为 $a^Ta$ 非可逆矩阵。



## 投影矩阵和最小二乘

[16. MIT线性代数---投影矩阵和最小二乘](https://zhuanlan.zhihu.com/p/45351889)

## 正交矩阵和Gram-Schmidt正交化法

[17. MIT线性代数---正交矩阵和Gram-Schmidt正交化法](https://zhuanlan.zhihu.com/p/46154870)

非常重要！

其中在Gram-Schmidt正交化法部分中，$b$ 在 $A$ 上的投影应为 $Ax$ 而非 $xA$，链接中公式有误

注意 $A,B$ 实际上是正交的向量，所以应采用在线上投影的公式，而非求可逆矩阵的方式

$B=b-A\hat x=b-A\frac{A^Tb}{A^TA}$

$A^TB=A^Tb-A^TA\frac{A^Tb}{A^TA}=0$



标准正交矩阵的优越性

对于标准正交矩阵 $Q$，正交向量常用 $q$ 表示，
$$
\begin{equation}
q_i\cdot q_j = \left\{
            \begin{array}{lcl}
            0 &\text{if} &i\neq j\\
            1 &\text{if} &i=j
            \end{array}  
       \right.
\end {equation}
$$
通过这个性质，得到 $Q^TQ=I$ ，对于投影矩阵 $P$，若 $A$ 为标准正交矩阵，则有 $P=Q(Q^TQ)^{-1}Q^T=QQ^T$ ，对于最小二乘法公式，有 $Q^TQ\hat x=Q^Tb\rightarrow \hat x=Q^Tb$

拆解开来，有 $\hat x_i=q_i^Tb$

## 行列式及其性质

[18. MIT线性代数---行列式及其性质](https://zhuanlan.zhihu.com/p/46157434)

补充对于4的证明

交换相等的两行，根据2得到 $-|A|=|A|\rightarrow |A|=0$

9证明了可逆的充要条件是行列式非0，否则的话 $|A^{-1}|=\frac{1}{|A|}$ 无意义

行列式为零，则为奇异矩阵

10表明，对行的性质同样对列存在，如一列为0，则行列式为0，因为可以将矩阵转置，行列式值不变（10的证明非常精髓，建议多看几遍）

比较重要的总结：

- **行列式的行**具有线性性质，见3，而整个行列式仅仅具有乘法性，不具有加法性，如 $|AB|=|A||B|,|A+B|\neq|A|+|B|$

## 行列式公式和代数余子式

[19. MIT线性代数---行列式公式和代数余子式](https://zhuanlan.zhihu.com/p/46158366)

求解行列式

1. 化简成三角矩阵后求解
2. 利用行列式的行满足线性关系的性质，将其一次拆分，即拆分成 $n!$ 项，这一公式的一大应用就是筛出每行每列都存在非零值的行列式，有时有奇效
3. Cofactor formula：对某一行（列）的每个元素一次展开 $(-1)^{i+j}|A_{n-1}|$ 

注意有时可以先化简再展开，或者对列进行展开（行列式的转置值不变）

对于视频中最后一题，求解 $\begin{vmatrix}1&1&0&0\\1&1&1&0\\0&1&1&1\\0&0&1&1\end{vmatrix}$

1. 化为上三角矩阵求解

2. 展开求解 $原式=\begin{vmatrix}1&0&0&0\\0&1&0&0\\0&0&1&0\\0&0&0&1\end{vmatrix}+\begin{vmatrix}1&0&0&0\\0&0&1&0\\0&1&0&0\\0&0&0&1\end{vmatrix}+\begin{vmatrix}0&1&0&0\\1&0&0&0\\0&0&1&0\\0&0&0&1\end{vmatrix}=-1$

3. Cofactor formula：对第一行展开 $原式=(-1)^{1+1}\begin{vmatrix}1&1&0\\1&1&1\\0&1&1\end{vmatrix}+(-1)^{1+2}\begin{vmatrix}1&1&0\\0&1&1\\0&1&1\end{vmatrix}=(-1)^{1+1}\begin{vmatrix}1&1&0\\1&1&1\\0&1&1\end{vmatrix}+(-1)^{1+2}(-1)^{1+1}\begin{vmatrix}1&1\\1&1\end{vmatrix}=|A_3|-|A_2|=-1$

   事实上，确实能够推出 $|A_n|=|A_{n-1}|-|A_{n-2}|$ ，该行列式满足六个一循环

## 克莱默法则、逆矩阵、体积

[20. MIT线性代数---克莱默法则、逆矩阵、体积](https://zhuanlan.zhihu.com/p/46159505)

$C$ 为代数余子式组成的矩阵，$C^T$ 称为原矩阵的伴随矩阵。

满足公式：$A^{-1}=\frac{1}{|A|}C^T$，只需证明 $AC^T=|A|I$

对于对角线上的元素，其均为Cofactor formula的展开，满足 $\sum_{k=1}^na_{ij}c_{jk}=|A|$，即依次对第 $i$ 行展开，对于非对角线上的元素，以第一行为例，$a_{11}c_{n1}+a_{12}c_{n2}+\dots+a_{1n}c_{nn}$ 等价于行列式在第 $n$ 行展开，即如 $\begin{vmatrix}a_{11}&a_{12}&\dots&a_{nn}\\a_{21}&a_{22}&\dots&a_{2n}\\\vdots&\vdots&\ddots&\vdots\\a_{11}&a_{12}&\dots&a_{1n}\end{vmatrix}$ 由于首行与第 $n$ 行相同，行列式为0，等式得证

基于如上公式，**克莱默法则**有：$x=A^{-1}b=\frac{1}{|A|}C^Tb=\frac{B}{|A|}$ ，其中 $B$ 是 $n$ 个行列式，每个行列式是将 $b$ 带入第 $i$ 列所得，如 $b_1=\begin{vmatrix}b_1&a_{12}&\dots&a_{nn}\\b_2&a_{22}&\dots&a_{2n}\\\vdots&\vdots&\ddots&\vdots\\b_n&a_{n2}&\dots&a_{nn}\end{vmatrix}$

行列式代表着对应空间的体积，其**证明思路**是：证明体积也满足行列式的三条基本性质

详细内容见链接内容，注意三角形的面积的行列式计算

## 特征值和特征向量

[21. MIT线性代数---特征值和特征向量](https://zhuanlan.zhihu.com/p/46716369)

对于方阵 $A$，有 $Ax=\lambda x$，对此 $\lambda$  称为特征值，$x$ 称为特征向量，即特征向量 $x$ 与 $Ax$ 指向统一方向

对于投影矩阵，若 b 在投影平面上，则 $Pb=b$，即 $\lambda=1,x=b$；若 b 为平面的法向量，则 $Pb=0$，即 $\lambda=0$，综上，投影矩阵的特征向量为 0，1

对于置换矩阵 $P=\begin{bmatrix}0&1\\1&0\end{bmatrix}$ ，$\begin{bmatrix}0&1\\1&0\end{bmatrix}\begin{bmatrix}1\\1\end{bmatrix}=\begin{bmatrix}1\\1\end{bmatrix}$ ，$\begin{bmatrix}0&1\\1&0\end{bmatrix}\begin{bmatrix}1\\-1\end{bmatrix}=-\begin{bmatrix}1\\-1\end{bmatrix}$ ，即置换矩阵 $P$ 的特征值为1，-1，特征向量为 $\begin{bmatrix}1\\1\end{bmatrix},\begin{bmatrix}1\\-1\end{bmatrix}$ ，注意这两个特征向量互相垂直

对于旋转矩阵（[旋转矩阵](https://www.cnblogs.com/WangGuiHandsome/p/10094784.html)），$Q=\begin{bmatrix}cos90&-sin90\\sin90&cos90\end{bmatrix}=\begin{bmatrix}0&-1\\1&0\end{bmatrix}$ ，其既为标准正交矩阵，又满足 $Q^T=-Q$ 这一极度不对称的特征，此时其特征值为纯虚数，即矩阵越是不对称，其特征值越偏向虚数

**注意：**

对于 $Ax=\lambda x,By=\sigma y$，二者的特征值不满足线性规律（即 $(A+B)X\neq(\lambda+\sigma)x,ABx\neq\lambda \sigma x$），因为往往特征向量不同，除非 $A,B$ 之中一矩阵为单位矩阵。

而对于A的幂次方，以二次方为例 $A^2x-\lambda Ax=\lambda^2 x$，特征值也为平方，特征向量不变

**规律：** 

- 特征值的和等于对角线上元素的和 $\sum_{1}^n \lambda_i=\sum_{1}^n a_{ii}$ ，对角线上元素的和称为迹
- 特征值的积等于矩阵的行列式 $\prod_i^n \lambda_i=|A|$



特征值及特征向量的求解见链接，主要思路 $Ax=\lambda x\rightarrow (A-\lambda I)x=0$ ，此时 x 处于零空间中，若 x 非零则 $rank(A-\lambda I)<n$ ，即其不可逆、奇异 $|A-\lambda I|=0$ ，结合其上规律求解，解得 $\lambda$ 后带入求解 $x$ 

## 对角化以及A的幂

[如何理解矩阵特征值？](https://www.zhihu.com/question/21874816)

[为什么实对称矩阵一定能对角化？](https://www.zhihu.com/question/38801697)

假设矩阵A具有n个线性无关的特征向量（不是所有矩阵都具有n个线性无关的特征向量），令其组成矩阵S，因此矩阵S必定可逆，有如下变换，其中 $\Lambda$ 为**对角特征矩阵**：
$$
AS=A\begin{bmatrix}x_1&x_2&\dots&x_n\end{bmatrix}
=\begin{bmatrix}\lambda_1x_1&\lambda_2x_2&\dots&\lambda_3x_n\end{bmatrix}
=\begin{bmatrix}x_1&x_2&\dots&x_n\end{bmatrix}\begin{bmatrix}\lambda_1&0&\dots&0\\0&\lambda_2&\dots&0\\\vdots&\vdots&\ddots&\vdots\\0&0&\dots&0\end{bmatrix}
=S\Lambda
$$
因此，我们可以得到公式 $S^{-1}AS=\Lambda$，即对角化，对角化的前提是 S 可逆，即 n 个特征向量独立

对比上一节对 $A^2$ 的理解，对角化后，$A^2=S\Lambda S^{-1}S\Lambda S^{-1}=S\Lambda^2 S^{-1}$ ，因为 $\Lambda$ 为对角矩阵，且特征值平方，故对焦特征矩阵也平方，**这就是特征值的应用，用于理解矩阵的幂次**

定理1，$A^k\rightarrow0\text{ as }k\rightarrow\infin\text{ if all }|\lambda_i|<1$

**KeyPoint：**存在 n 个独立的特征向量 $\Longleftrightarrow$ 所有特征值不同，事实上，在有特征值相同的情况下，也可能存在 n 个独立的特征向量，例如单位矩阵特征值都为1，但是其存在 n 个独立的特征向量，也存在特征向量不足 n 个的情况 （[如何理解不同特征值对应的特征向量线性无关](https://zhuanlan.zhihu.com/p/30454490)）

对于 $u_{k+1}=Au_k$ ，有 $u_{k}=A^ku_0$，以斐波那契数列为例，$F_0=0,F_1=1,F_{k+2}=F_{k+1}+F_k$，解法重点，添加一行 $F_{k+1}=F_{k+1}$，得到方程组 $\begin{equation}
\left\{
\begin{aligned}
F_{k+2}&=F_{k+1}+F_k\\
F_{k+1}&=F_{k+1}\\
\end{aligned}
\right.
\end{equation}$ ，令 $u_k=\begin{bmatrix}F_{k+1}\\F_k\end{bmatrix}$，得到如下形式
$$
\begin{bmatrix}F_{k+2}\\F_{k+1}\end{bmatrix}=\begin{bmatrix}1&1\\1&0\end{bmatrix}\begin{bmatrix}F_{k+1}\\F_{k}\end{bmatrix}\Longrightarrow u_{k+1}=\begin{bmatrix}1&1\\1&0\end{bmatrix}u_k
$$
对 $A=\begin{bmatrix}1&1\\1&0\end{bmatrix}=S\Lambda S^{-1}$ ，A的特征值为 $\lambda_1=\frac{1+\sqrt{5}}{2},\lambda_2=\frac{1-\sqrt{5}}{2}$，特征向量为 $\begin{bmatrix}\lambda\\1\end{bmatrix}\Rightarrow x_1=\begin{bmatrix}\lambda_1\\1\end{bmatrix},x_2=\begin{bmatrix}\lambda_2\\1\end{bmatrix}$，$S=[x_1,x_2]$，可以发现这里的特征值其实就是高中数列所求的特征值，由于 $x_1,x_2$ 线性无关，可以用**来表示任意二维向量**，因此 $F_0=\begin{bmatrix}1\\0\end{bmatrix}=c_1x_1+c_2x_2=Sc$，即
$$
u_{k+1}=\begin{bmatrix}1&1\\1&0\end{bmatrix}u_k=A^ku_0=S\Lambda^kS^{-1}Sc=\Lambda^kSc\\或\quad u_{k+1}=\begin{bmatrix}1&1\\1&0\end{bmatrix}u_k=A^ku_0=\lambda_1^kx_1c_1+\dots+\lambda_n^kx_nc_n=\Lambda^kSc
$$
根据定理1，$|\lambda_2|<1$，则增长速度主要由 $\lambda_1$ 决定，即 $u_k\approx\lambda^k x_1c_1$，换句话说，随着 k 的增大，**矩阵越来越贴近最大的特征值**

## 微分方程和exp(At)

[线性代数之——微分方程和 exp(At)](https://zhuanlan.zhihu.com/p/51128280)

## 马尔科夫矩阵；傅里叶级数

马尔科夫矩阵的条件：

1. 所有元素大于等于零，因为马尔科夫矩阵与概率相关，概率非负
2. 每列的和为1

示例：$A=\begin{bmatrix}0.1&0.01&0.3\\0.2&0.99&0.3\\-.7&0&0.6\end{bmatrix}$ 

马尔科夫矩阵的特点：

1. 特征值中一定包含1，确保其具有稳态
   $$
   A-1I=\begin{bmatrix}-0.9&0.01&0.3\\0.2&-0.01&0.3\\-.7&0&-0.6\end{bmatrix}
   $$
   证明其行列式为0有两种方法

   - 由于马尔科夫矩阵每列和为1，所以 $A-1I$ 每列和为0，同时加到最后一行得到全零行，此时行列式为0
   - 同样基于其列和为0，我们可以得到向量 $\begin{bmatrix}1\\1\\1\end{bmatrix}$ 存在其左零空间中，因此该矩阵行线性相关，行列式为0

   **BY THE WAY** ：$A$ 的特征向量等同于 $A^T$ 的特征向量，因为 $|A|=|A^T|$ ，所以 $|A-\lambda I|=|(A-\lambda I)^T|=|A^T-\lambda I|$，即二者特征值相等

2. 其他所有特征值绝对值小于1，$|\lambda_i|<1$，确保其他项收敛

马尔科夫矩阵的应用：稳态模型，稳态最终有特征为1的部分决定



标准正交基的投影问题

对于一组标准正交基 $q_1,q_2,\dots,q_n$ ，有 $v=x_1q_1+x_2q_2+\dots+x_nq_n\Longrightarrow v=Qx$ ，求解 $x_1,x_2,\dots,x_n$

只需依次左乘正交基即可：$q_1^Tv=x_1q_1^Tq_1+x_2q_1^Tq_2+\dots=x_1q_1^Tq_1=x_1\Longrightarrow x=Q^{-1}v=Q^Tv$ 

在傅里叶级数的应用：

$f(x)=a_0+a_1cosx+b_1sinx+a_2cos2x+b_2sin2x+\cdots$

函数的正交性，对于函数 $f$ 以及 $g$，$f^Tg=\int f(x)g(x)dx$，由于 $f(x)=f(x+2\pi)$ ，对于其中涉及的三角函数，两两正交，傅里叶级数中，使用正交函数来替代正交向量

因此 $\int_0^{2\pi}f(x)cosxdx=a_1\int_0^{2\pi} cos^2xdx=\pi\Longrightarrow a_1=\frac{1}{\pi}\int_0^{2\pi}f(x)cosxdx$ ，其他项的求解以此类推

## 复习二



## 对称矩阵及其正定性

- 实对称矩阵的特征值一定是实数（[实对称矩阵的特征值为实数的Proof](http://blog.sina.com.cn/s/blog_4cb6ee6c0100g56l.html)），如果是复数矩阵，需要满足 $\bar A^T=A$，上述结论才成立（可从证明过程中理解）

- 对称矩阵的特征向量一定是正交的（[对称矩阵的特征向量两两正交的证明](https://blog.csdn.net/bingecuilab/article/details/47209037)），因此在对称矩阵中 $A=S\Lambda S^{-1}=Q\Lambda Q^{-1}=Q\Lambda Q^T$

$A=Q\Lambda Q^T=\lambda_1q_1q_1^T+\lambda_2q_2q_2^T+\dots+\lambda_nq_nq_n^T $ 其中 $q_iq_i^T$ 为投影矩阵，因此**每个对称矩阵都可以表示成若干正交的投影矩阵的线性组合**



**数值计算的小定理**：主元符号的分布同特征值符号的分布相同，即主元中10个为正，10个为负，那么特征值中也是10个为正，10个为负（这在微分方程中很有用，因为只需知道符号就能知道稳定与否）

**正定矩阵（positive definite symmetric matrix）**：不完整，见如下笔记

- 所有特征值为正数 $\Longleftrightarrow$ 所有主元为正 $\Longleftrightarrow$ 行列式为正
- 所有子行列式为正

## 复数矩阵和快速傅里叶变换



## 正定矩阵和最小值

正定矩阵的判断条件：

- $\lambda_1>0,\quad \lambda_2>0$
- $a>0,\quad ac-b^2>0$
- pivots $a>0,\quad\frac{ac-b^2}{a}>0$
- :star: $x^TAx>0$

**半正定矩阵**：对于任意非零实向量 x，都有 $xAx^T\geq 0$

**从代数角度理解正定矩阵**：

以二阶矩阵为例，对于矩阵 $\begin{bmatrix}2&6\\6&20\end{bmatrix}$ 为例，

$$
\begin{bmatrix}x_1\\x_2\end{bmatrix}\begin{bmatrix}2&6\\6&20\end{bmatrix}\begin{bmatrix}x_1&x_2\end{bmatrix}=2x_1^2+12x_1x_2+20x_2^2
\\ 
原式=
\begin{bmatrix}x_1\\x_2\end{bmatrix}\begin{bmatrix}2&6\\0&2\end{bmatrix}\begin{bmatrix}x_1&x_2\end{bmatrix}=
2(x_1+3x_2)^2+2x_2^2
$$
注意化简后的矩阵即为化简后的代数表达式，其中2，2为主元，通过化简后（配方后）的表达式，我们能够判断是否是正定矩阵（是否恒大于0）

**判断极小值**：$f(x_1,x_2,\dots,x_n)$ 存在极小值的条件是二阶导数矩阵，如 $\begin{bmatrix}f_{xx}&f{xy}\\f_{yx}&f_{yy}\end{bmatrix}$ 是正定矩阵，这也就是微积分中求极值时要求 $f_{xx}f_{yy}-f_{xy}f_{yx}>0$ (即二阶导数矩阵的行列式>0)的原因，即要求其为正定矩阵

**从图形角度理解正定矩阵**：

对于矩阵 $A=\begin{bmatrix}2&-1&0\\-1&2&-1\\0&-1&2\end{bmatrix}$ ，A的子矩阵的行列式分别为 $2,3,\dots,n$ ，通过特征值的积等于行列式的值，$\lambda=2,\frac{3}{2},\dots,\frac{n+1}{n}$ ，其特征值同时也是其主元，对于A而言，其特征向量为 $2-\sqrt{2},2,2+\sqrt{2}$ 
$$
f=xAx^T=2x_1^2+2x_2^2+2x_3^2-2x_1x_2-2x_2x_3
$$
如果取 $f=c$ 将会得到一个变形椭球体，其中三个特征向量分别指向长中短轴的方向，特征值的大小决定轴的长度 

[主轴定理 (Principal axis theorem)](https://blog.csdn.net/linkequa/article/details/88757280)

## 相似矩阵和若尔当形

对于矩阵 $A_{m\times n}$，$A^TA$ 一定是正定矩阵，因为 $x^TA^TAx=(Ax)^TAx=\Vert Ax\Vert>0$，rank=n

因此如果 $A^TA$ 可逆，那么最小二乘法就可行

**矩阵相似**：若存在 M，使得 $B=M^{-1}AM$ ，则称矩阵 A、B 相似，如 $S^{-1}AS=\Lambda$ 则称 A 相似于 $\Lambda$ 

两矩阵相似，他们具有相同的特征值
$$
Ax=\lambda x\Longrightarrow AMM^{-1}x=\lambda x\Longrightarrow 
M^{-1}AMM^{-1}x=\lambda M^{-1}x\Longrightarrow BM^{-1}x=\lambda M^{-1}x
$$
由此可知 A、B 特征值相同，$M^{-1}x$ 为 B 的特征向量

 **若尔当定理**：任意矩阵 A 都相似于一个 若尔当矩阵 J, $J=\begin{bmatrix}J_1&&&\\&J_2&&\\&&\ddots&\\&&&J_d\end{bmatrix}$ ，$J_i$ 为若尔当块

若尔当块形如 $J_i=\begin{bmatrix}\lambda_1&1&&\\&\lambda_2&1&\\&&\ddots&\\&&&\lambda_d\end{bmatrix}$ ，每个块中只包含一个特征值，即（`len(unique(lambda))=1`）因此一个矩阵若尔当块的数目等于特征值的数目，若尔当块不同，矩阵不相似
$$
A=\begin{bmatrix}0&1&0&0\\0&0&1&0\\0&0&0&0\\ 0&0&0&0\end{bmatrix}
=\begin{bmatrix}J1_{3\times3}&\\&J2_{1\times1}\end{bmatrix}\qquad 
B=\begin{bmatrix}0&1&0&0\\0&0&0&0\\0&0&0&1\\0&0&0&0\end{bmatrix}
=\begin{bmatrix}J1_{2\times2}&\\&J2_{2\times2}\end{bmatrix}
$$
如图由于若尔当块大小不同，所以矩阵A、B不相似

## 奇异值分解

[奇异值的物理意义](https://www.zhihu.com/question/22237507)

[SVD和PCA](https://www.zhihu.com/question/38319536)

SVD singular value decomposition 即对于任意矩阵A，有 $A=U\Sigma V^T$ ，其中 $U、V$ 为正定矩阵，$\Sigma$ 为对角矩阵



## 线性变换及对应矩阵



## 基变换和图像压缩



## 复习三



## 左右逆和伪逆



## 期末复习

[互逆矩阵的特征向量相同，特征值为倒数](https://zhidao.baidu.com/question/648581912836061085.html)