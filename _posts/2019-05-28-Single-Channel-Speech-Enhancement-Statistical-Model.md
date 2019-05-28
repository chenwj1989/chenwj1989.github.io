---
layout: post
title:  "单通道语音增强之统计信号模型"
date:   2019-05-28 14:32:04 +0700
categories: [speech, 语音]
---

<br/> 
*  目录
{:toc}

## 1. 信号估计理论简述

信号估计理论是现代统计处理的基础课题[@ZhangXianDa2002ModernSP]，在通信、语音、图像领域均有广泛应用。语音增强，就是从带噪的语音测量信号中估计原始的无噪语音，这是典型的信号估计问题。
《语音增强–理论与实践》[@loizou2007speech]一书中列举了用于语音增强的一系列统计模型。

假设麦克风采集到的带噪语音序列为$$y[n]$$，并且噪声都是加性噪声。则带噪语音序列为无噪语音序列与噪声序列的和。原始语音信号与噪声均可视为随机信号。

 $$y[n] = x[n] + d[n]$$

在时域对$x[n]$进行估计是非常困难的，通过傅立叶变换，我们可以将信号分解为频域上互相独立的系数。信号估计模型转变为对每一个频点的系数进行估计的模型，不同频点之间的参数是相互独立的。

$$Y(\omega_{k}) = X(\omega_{k}) + D(\omega_{k})$$

这个方法就叫做统计信号谱分析（Statistical Spectral
Analysis)。显然地，纯净信号谱$X(\omega_{k})$带有幅度与相位两参数，我们实际上是对幅度$X_{k}$和相位$\theta_{y}(k)$进行参数估计。

$$Y_{k}e^{j\theta_{y}(k)} = X_{k}e^{j\theta_{x}(k)} + D_{k}e^{j\theta_{d}(k)}$$

重组谱幅度和谱相位估计值即可恢复纯净语音谱，估计值用上标来表示：$\hat{X}(\omega_{k})=\hat{X}(k)e^{j\hat{\theta}_{x}(k)}$。

实际信号的幅度和相位是不方便直接用在运算过程中的，因为信号取值范围不定，且瞬时变化。在噪声抑制领域，更常用的语音谱估计方法是对抑制增益(Suppression
Gain)进行估计，不同的估计准则称为抑制准则(Supression Rule)。

$$\hat{X}(\omega_{k}) = H_{k}Y(\omega_{k})$$

通常会根据先验信噪比、后验信噪比来估计抑制增益$H_{k}$。
并且可以在只有噪声出现的时刻更新$H_{k}$，
在语音存在的时刻进行抑制，无须每帧去调用噪声抑制算法，计算过程比直接估计信号谱灵活。

综上，语音增强的典型流程就是：

1.  对带噪语音y\[n\]分帧， 每一帧进行DFT得到$Y(\omega_{k})$。

2.  估计或者沿用上一帧的抑制增益$H_{k}$，得到纯净语音谱$\hat{X}(\omega_{k})$

3.  对$\hat{X}(\omega_{k})$进行IDFT,得到纯净语音序列的估计$x[n]$。

为了估计模型的建模，对测量信号、估计信号、噪声信号都需要作一些数学上的假设和简化。其中对噪声一般会作以下假设：

1.  噪声是与语音独立的加性噪声；

2.  每一帧噪声的统计分布是稳态的；

3.  噪声的傅立叶级数是零均值复高斯分布。

$$\left[D(\omega_{k})\right] = 
\frac{1}{\pi\lambda_{d}(k)}
\exp\left[-\frac{|D(\omega_{k})|^2}{\lambda_{d}(k)}\right]$$

## 2. 最大似然估计ML

基如果不考虑信号的先验分布，即认为信号值是确定信号，而不是随机信号，我们只需要分析含有信号$x$为参数的带噪信号$y$的概率分布$p(y;x)$，并使之最大。这种估计方法叫最大似然估计器（Maximum
Likelihood
Estimation）。将$y$的带参概率分布$p(y;x)$称为似然函数（Likelihood
Function）。对纯净信号$x$的估计，表达为求解合适的$x$值，使得似然函数$p(y;x)$最大。
$$\hat{x} = \mathop{\arg\max}_{x}p(y;x)$$

文献[@mcaulay1980speech]最早将最大似然估计法用在语音增强领域。对于纯净语音，可以假设纯净语音幅度$X_{k}$和相位$\theta_{y}(k)$是未知但确定，无需考虑其先验概率分布。最大似然语音增强模型表达为：

$$\hat{X}_{k},\hat{\theta}_{k} = \mathop{\arg\max}_{X_{k},\theta_{k}} p\left[Y(\omega_{k});X_{k},\theta_{k}\right]$$

把$D_{k}e^{\theta_{d}(k)}=Y_{k}e^{\theta_{y}(k)} - X_{k}e^{\theta_{x}(k)} $代入噪声零均值复高斯分布公式中，得到：

$$p\left[Y(\omega_{k});X_{k},\theta_{k}\right] = \frac{1}{\pi \lambda_{d}(k)}\exp\left[-\frac{|Y(\omega_{k}) - X_{k}e^{j\theta_{x}(k)}|^2}{\lambda_{d}(k)}\right]$$

$$H_{k} = \frac{1}{2}+\frac{1}{2}\sqrt{\frac{\xi_{k}}{1+\xi_{k}}}$$

## 3.贝叶斯估计

如果比最大似然估计更进一步，考虑待估计量$x$也是随机变量，且$x$的先验分布为$p(x)$，这种假设下的估计方法叫做贝叶斯估计[@ZhangXianDa2002ModernSP]。定义估计值$\hat{x}$与实际值$x$之间的误差函数为$c(\hat{x},x)$，贝叶斯估计器的目标即为找出是平均误差$E[c(\hat{x},x)]$最小的估计值$x$。

$$\hat{x} = \mathop{\arg\min}_{\hat{x}} E[c(\hat{x}, x)]$$

对于待估计的纯净语音谱，贝叶斯估计器可以表达为：

$$\hat{X}(\omega_{k}) = \mathop{\arg\min}_{\hat{X_{k}}} E[c(\hat{X}(\omega_{k}), X(\omega_{k}))]$$

误差$c(\hat{x},x)$的期望值取决于待测信号与测量信号的联合概率分布。

$$E[c(\hat{x}, x)] = \int_{x}\int_{y}c(\hat{x}, x)p(x, y)dxdy$$

对联合概率密度进行分解，提取$p(y)$与条件概率密度$p(x|y)$：

$$E[c(\hat{x}, x)] = \int_{y} \left[ \int_{x}c(\hat{x}, x)p(x|y)dx\right] p(y)dy$$

因为估计值$\hat{x}$与$p(y)$相互独立，对$\hat{x}$的求解可以等价于对$\int_{x}c(\hat{x}, x)p(x|y)dx$求极大值，也就是：

$$\hat{x} = \mathop{\arg\min}_{\hat{x}} E[c(\hat{x}, x)|y]$$

误差函数$c(\hat{x}, x)$的计算模型，会引出不同种类的估计器。典型的误差函数有几种类型[@ZhangXianDa2002ModernSP]：

1.  平方误差函数，对应最小均方估计。

    $$c(\hat{x}, x)  = (\hat{x} - x)^2$$

2.  绝对值误差函数，对应条件中位数估计。

    $$c(\hat{x}, x)  = |\hat{x} - x|$$

3.  均匀误差函数，对应最大后验估计。

 $$\begin{equation}c(\hat{x}, x)   =
            \left\{
            \begin{array}{lr}
            1 \quad (|\hat{x}-x|\geq\Delta/2) &  \\
            &\Delta>0\\
            0 \quad (|\hat{x}-x|<\Delta/2) & 
            \end{array}
            \right.\end{equation}$$

### 3.1 最小均方估计（MMSE）

最小均方估计(Minimum Mean Square Error
Estimation)[@ephraim1985speech]使用平方误差函数，平均误差为：

$$E[c(\hat{x}, x)|y]=\int_{x}(\hat{x} - x)^2 p(x|y)dx$$

为求平均误差的极大值，可对估计量$\hat{x}$求导，并求极值点。

$$\frac{\partial E[c(\hat{x}, x)|y]}{\partial \hat{x}}=\int_{x}2(\hat{x} - x) p(x|y)dx =0$$

显然极值点的$\hat{x}$为被估计量$x$的条件均值：

$$\hat{x} =\int_{x}x p(x|y)dx = E[ x| y]$$

亦可以使用贝叶斯公式展开，得到： 

$$\begin{split}
    \hat{x} = \int_{x}x \frac{p(y|x)p(x)}{p(y)} dx \\
            = \frac{\int_{x}xp(y|x)p(x)dx}{\int_{x}p(y|x)p(x)dx}
    \end{split}$$

在语音增强的模型里，纯净语音谱的估计为其在带噪语音谱下的条件均值。

$$\hat{X}(\omega_{k}) = E[X(\omega_{k})| Y(\omega_{k})]$$

上式中，为得到纯净语音谱需要分别估计谱幅度$\hat{X}_{k}$和谱相位$\hat{\theta}_{x}(k)$。

$$\hat{X}_{k}e^{j\hat{\theta}_{x}(k)} = E[ X_{k}e^{j\theta_{x}(k)}| Y(\omega_{k})]$$

同时估计谱幅度和谱相位是很难的，研究者提出了许多分别估计谱幅度和谱相位的方法，估计完成后再用两者重组复语音信号。

$$\hat{X}_{k} = E[ X_{k}| Y(\omega_{k})]$$

$$e^{j\hat{\theta}_{x}(k)} = E[ e^{j\theta_{x}(k)}| Y(\omega_{k})]$$

#### 3.1.1 MMSE谱幅度估计

最小均方根估计器MMSE short-time spectral amplitude

$$\hat{X_{k}} = E[ X_{k}| Y(\omega_{k})]$$

$$v_{k} = \frac{\xi_{k}}{1+\xi_{k}}\gamma_{k}$$

$$H_{k} =\frac{\sqrt{\pi v_{k}}}{ 2\gamma_{k}}[(1+v_{k})I_{0}(\frac{v_k}{2})+v_{1}I_{1}(\frac{v_{k}}{2})]\exp(\frac{-v_{k}}{2})$$

```python
	def mmse_stsa_gain(parameters=None):
		gamma = parameters['gamma']
		ksi = parameters['ksi']
			
		vk = ksi * gamma / (1 + ksi)
		j0 = np.i0(vk / 2)
		j1 = np.i1(vk / 2)
			
		A = sqrt(pi * vk) / 2 / gamma
		B = (1 + vk) * j0 + vk * j1
		C = np.exp(-0.5 * vk)
		gain = A * B * C
		return gain
```

#### 3.1.2 MMSE对数谱幅度估计

对数最小均方根估计器The MMSE log spectral amplitude (MMSE-LSA)，
或者缩写为LogMMSE估计器。

$$c(\hat{X_{k}}, X_{k})  = (\log{X_{k}}- \log{X_{k}})^2$$

$$\log{\hat{X_{k}}} = E[ \log{X_{k}}| Y(\omega_{k})]$$

$$H_{k} = \frac{\xi_{k}}{1+\xi_{k}}\exp(\frac{1}{2}\int_{v_{k}}^{\infty}\frac{e^{-t}}{t}dt)$$

```python
	def logmmse_gain(parameters=None):
		gamma = parameters['gamma']
		ksi = parameters['ksi']
		A = ksi / (1 + ksi)
		vk = A * gamma
		ei_vk = 0.5 * expn(1, vk)
		gain = A * np.exp(ei_vk)
		return gain
```

#### 3.1.3 MMSE平方谱幅度估计

频谱幅度平方估计器MMSE magnitude squared[@wolfe2003efficient]

$$\hat{X_{k}^2} = E[ X_{k}^2| Y_{k}]$$

$$\hat{X_{k}^2} = \frac{\xi_{k}}{1+\xi_{k}} (\frac{1+v_{k}}{\gamma_{k}})Y_{k}^2$$

$$H_{k} = \sqrt{\frac{\xi_{k}}{1+\xi_{k}} (\frac{1+v_{k}}{\gamma_{k}})}$$

```python
	def mmse_sqr_gain(parameters=None):
		gamma = parameters['gamma']
		ksi = parameters['ksi']
			
		vk = ksi * gamma / (1 + ksi)
		j0 = np.i0(vk / 2)
		j1 = np.i1(vk / 2)
			
		A = ksi / (1 + ksi)
		B = (1 + vk) / gamma
		gain = sqrt(A * B)
		return gain
```

### 3.2 最大后验估计 MAP

当贝叶斯估计器采用均匀误差函数时，平均误差为： 

$$\begin{split}
    E[c(\hat{x}, x)|y] &= \int_{-\infty}^{\hat{x}-\Delta/2}  p(x|y)dx + \int_{\hat{x}+\Delta/2}^{\infty}  p(x|y)dx   \\
                       &= 1-\int_{\hat{x}\Delta/2}^{\hat{x}+\Delta/2}  p(x|y)dx 
\end{split}$$

显然要使得平均误差最小，就是要求目标估计$\hat{x}$，使得
$p(x|y)$最大。这种估计模型称作最大后验估计(Maximum A Posteriori
Estimation)。这个模型的意思是只有估计值$\hat{x}$与原始值$x$相等，误差才为0，其他时候误差均匀为1。估计值可以表达为：

$$\hat{x} =\mathop{\arg\min}_{x} p(x|y)$$

在语音增强的模型里，纯净语音谱的估计为其在带噪语音谱下的条件均值。

$$\hat{X}(\omega_{k}) = \mathop{\arg\min}_{ X(\omega_{k})}\quad p\left[X(\omega_{k})| Y(\omega_{k})\right]$$

文献[@wolfe2003efficient]提出了两种基于最大后验估计(Maximum A Posteriori
Estimation)的语音增强算法。一种是同时求解幅度和相位的混合最大后验估计：
另外一种是单纯估计幅度的方法，两种估计器最后的噪声抑制增益略有不同。

#### 3.2.1 幅度和相位混合最大后验估计

[@wolfe2003efficient]

$$H_{k} = \frac{\xi_{k}+\sqrt{\xi_{k}^2+2(1+\xi_{k})\frac{\xi_{k}}{\gamma_{k}}}}{2(1+\xi_{k})}$$

```python
def map_joint_gain(parameters=None):
	gamma = parameters['gamma']
	ksi = parameters['ksi']
	
	eps = 1e-6
	gain = (ksi + sqrt(ksi^ 2 + 2 * (1.0 + ksi)* ksi/ (gamma + eps))) / 2.0/ (1.0 + ksi)
	return gain
```

#### 3.2.2 纯幅度最大后验估计

[@wolfe2003efficient]

$$H_{k} = \frac{\xi_{k}+\sqrt{\xi_{k}^2+(1+\xi_{k})\frac{\xi_{k}}{\gamma_{k}}}}{2(1+\xi_{k})}$$

```python
def map_sa_gain(parameters=None):
	gamma = parameters['gamma']
	ksi = parameters['ksi']
	
	eps = 1e-6
	gain = (ksi + sqrt(ksi^ 2 + (1.0 + ksi)* ksi/ (gamma + eps))) / 2.0/ (1.0 + ksi)
	return gain
```
<br/>




























































