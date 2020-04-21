# MIT线性代数

课程：[麻省理工线性代数](https://www.bilibili.com/video/BV1Kb411M72F)

## 方程组的几何解释

从行的角度看，矩阵可以看作是方程组（方程组也可以看成面，于是化为面的交点）；

从列的角度，每一列可以看作一个三维向量，化为三个向量的线性和，于是对于如下问题

:star: `Can I solve Ax = b for every b ?` $\Longleftrightarrow$ `Do the linear combinations of the columns fill the 3-D space?`

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

