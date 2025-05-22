# Crystal plasticity constitutive model framework

# Crystal deformation

Deformation gradient tensor $\mathbf{F}$ and velocity gradient $\bf{l}$ ：

$$
\bf{F}=\frac{\partial\bf{{x}}}{\partial\bf{{X}}}
$$

$$
\bf{l}=\frac{\partial \bf{v}}{\partial \bf{X}}
$$

The relationship between them is：

$$ \dot{\mathbf{F}}=\mathbf{l}\mathbf{F} $$

$$ \mathbf{F}=\mathbf{F}^e\mathbf{F}^p $$

$$ \mathbf{l}=\dot{\mathbf{F}}\mathbf{F}^{-1}=\dot{\mathbf{F}^e}(\mathbf{F}^e)^{-1}+\mathbf{F}^e\dot{\mathbf{F}^p}(\mathbf{F}^p)^{-1}(\mathbf{F}^e)^{-1}$$

$$ \mathbf{l}^e=\dot{\mathbf{F}^e}(\mathbf{F}^e)^{-1}$$

$$ \mathbf{l}^p=\mathbf{F}^e\dot{\mathbf{F}^p}(\mathbf{F}^p)^{-1}(\mathbf{F}^e)^{-1}$$

For every velocity gradient, there exists a symmetric part (strain rate tensor) and an antisymmetric part (spin tensor)：

$$ \mathbf{l}=\mathbf{d}+\mathbf{w} $$

$$ \mathbf{l}^e=\mathbf{d}^e+\mathbf{w}^e $$

$$ \mathbf{l}^p=\mathbf{d}^p+\mathbf{w}^p $$

At each solution step, the deformation-related tensors of the grain (strain, orientation matrix) are updated.：

$$ \mathbf{\varepsilon}_{t+\Delta t}=\mathbf{\varepsilon}_t+\mathbf{d}\Delta t $$

$$ \mathbf{M}\_{t+\Delta t}=\mathbf{M}\_{t}(\mathbf{I}+\mathbf{w}^e\Delta t)^T $$

Here, the orientation matrix of the grain is:

$$ \mathbf{M}=\begin{bmatrix} u & m & h \\
v & n & k \\
w & o & l \end{bmatrix} $$

The x, y, and z axes of the reference coordinate system correspond respectively to the following directions in the grain coordinate system $(uvw),(mno),(hkl)$

For the elastic deformation gradient, its polar decomposition yields a rotation matrix and a stretch (deformation) matrix.：

$$
\mathbf{F}^e=\mathbf{RU}
$$

则有：

$$
\mathbf{M}=\mathbf{M}_0\mathbf{R}^T
$$

# Elastoplastic constitutive model

$$
\dot{\mathbf{\sigma}}-\mathbf{w}^e\mathbf{\sigma}+\mathbf{\sigma}\mathbf{w}^e+\mathbf{\sigma} tr(\mathbf{d}^e)=\mathbb{C}:\mathbf{d}^e
$$

$$ \dot\sigma=\mathbf{C}(d-d_p)-\sigma tr(d)+w_e\sigma-\sigma w_e $$

$$ =\mathbf{C}d-\sigma tr(d) + w\sigma - \sigma w - \mathbf{C}d_p - w_p\sigma + \sigma w_p $$

### Convergent solution: Newton–Raphson iteration method

$$
\begin{align*}f(X)&=(\mathbf{C}d-\sigma tr(d) + \mathbf{\Sigma} w - \mathbf{C}d_p -   \mathbf{\Sigma} w_p - \dot{\sigma})dt \\ 
&=(\mathbf{C}d-\sigma tr(d) + \mathbf{\Sigma} w - \mathbf{C}d_p -   \mathbf{\Sigma} w_p)dt - \Delta\sigma\end{align*} 
$$

$$ X=(w,d,\Delta\sigma) $$

Where：

$$
\mathbf{\Sigma}=\begin{pmatrix} 0 & \sigma_5 & \sigma_6 \\\\
\sigma_4 & 0 & -\sigma_6 \\\\
-\sigma_4 & -\sigma_5 & 0 \\\\ 
\frac12(\sigma_3-\sigma_2) & -\frac12\sigma_6 & -\frac12\sigma_5\\\\
-\frac12\sigma_6 & \frac12(\sigma_3-\sigma_1) & \frac12\sigma_4\\\\
\frac12\sigma_5 & \frac12\sigma_4 & \frac12(\sigma_2-\sigma_1) \end{pmatrix}
$$

Thus, the gradient becomes：

$$
\frac{\partial f(X)}{\partial X} = \mathbf{C'} \partial d+\mathbf{\Sigma}\partial w +\frac{\partial f(X)}{\partial{\Delta\sigma}}\partial{\Delta\sigma}
$$

$$
\mathbf{C'}=\mathbf{C}-\begin{pmatrix}\sigma&\sigma&\sigma&0&0&0\end{pmatrix}
$$

To be consistent with the code, the tensor is written in 6-dimensional form.：

$$ \frac{\partial f(X)}{\partial \Delta \sigma}=-I-C\frac{\partial\tilde{d_p}}{\partial \Delta\sigma}-\Sigma N^{-1}\frac{\partial\tilde{w_p^*}}{\partial \Delta\sigma} $$

$$ =-I-C(M_{cpl}\frac{\partial\tilde{d_p}^c}{\partial \Delta\sigma^c}M_{cpl}^T)-\Sigma N^{-1}(M_{cpl}\frac{\partial\tilde{d_p^{*c}}}{\partial \Delta\sigma^c}M_{cpl}^T)_{4:6} $$

Where：

$$ \tilde{d_p} = \begin{bmatrix} d_1 & d_2 &d_3 & 2d_{23} & 2d_{13} & 2d_{12} \end{bmatrix}] $$

$$ \tilde{w_p^*} = \begin{bmatrix} 0 & 0 & 0 & 2w_{23} & 2w_{13} & 2w_{12} \end{bmatrix}] $$

The superscript c denotes the crystal coordinate system, and $M_{cpl}$ is the corresponding 6×6 rotation matrix for the compliance tensor。

$$ S = \begin{pmatrix} s_1 n_1 & s_2 n_2 & s_3 n_3 & (s_2 n_3 + s_3 n_2) & (s_1 n_3 + s_3 n_1) & (s_1 n_2 + s_2 n_1) \end{pmatrix}^T $$

$$ A = \begin{pmatrix} 0 & 0 & 0 & (s_2 n_3 - s_3 n_2) & (s_1 n_3 - s_3 n_1) & (s_1 n_2 - s_2 n_1) \end{pmatrix}^T $$

$$ \frac{\partial{d_p^c}}{\partial \Delta\sigma^c}=SS^T\frac{\partial\dot\gamma dt}{\partial\tau} $$

$$ \frac{\partial{w_p^c}}{\partial \Delta\sigma^c}=AS^T \frac{\partial\dot\gamma dt}{\partial\tau} $$

Depending on the boundary condition settings, there are

$$
\begin{pmatrix} d \\\\ w \\\\ \Delta \sigma \end{pmatrix}\_{15\times 1} 
= \begin{pmatrix} I\_{3\times 3} & 0 & 0 & 0\_{3\times 6} \\\\ 0 & 0.5I\_{3\times 3} & 0.5I\_{3\times 3} & 0\_{3\times 6} \\\\
0 & 0.5I\_{3\times 3} & -0.5I\_{3\times 3} & 0\_{3\times 6} \\\\
0 & 0 & 0 & I\_{6\times 6} \end{pmatrix}\_{15\times 15} \begin{pmatrix} L \\\\ 
\Delta \sigma \end{pmatrix}\_{15*1}
$$

Derive the objective function $f(X)=0$ and the gradient $\partial f(X)/\partial X$ After that, $X$ can be solved using the Newton-Raphson method, that is

$$
X_{n+1}= X_n + \frac{f(X)}{\partial f(X)/\partial X}
$$

### When Newton-Raphson fails to converge: apply a Newton descent (line search) method.

$$
X_{n+1}= X_n + c\cdot\frac{f(X)}{\partial f(X)/\partial X}
$$

- If $|f_{n+1}|/|f_n|>1$ ，Let $c=1/2c$ ，and recompute the current step；if $|Y_{n+1}|/|Y_n|\le1$ ，Let $c=1$ ；

# Plastic deformation mechanis

## Slip

Computation of the plastic velocity gradient

$$ \mathbf{l}^p=\Sigma_\alpha\dot\gamma^\alpha \mathbf{s}^\alpha \mathbf{n}^{\alpha T}  $$

$$ \mathbf{s}^\alpha=\mathbf{F}^e\mathbf{s}_0^\alpha $$

$$ \mathbf{n}^\alpha=(\mathbf{F}^{e-1})^{T}\mathbf{n}_0^\alpha $$

Dislocation velocity and resolved shear stress

$$
v=L/(t_w + t_r)
$$

At present, there is no universally accepted form for the waiting time $t_w$ and running time $t_r$ in the expression. Here, we consider the effects of two types of forces: strong pinning forces and weak pinning forces. The strong pinning force is temperature-independent, while the weak pinning force is temperature-dependent. Thus

$$
\begin{align*} & t_w= \frac{1}{\nu_0} \exp(\frac{Q_a}{k_BT}), \\\\
&Q_a = Q_0[1-\text{sgn}(|\tau|-\tau_f)(\frac{||\tau|-\tau_f|}{\tau_c})^\xi] \\\\ 
&t_r=L \left(v_{s}\left[\sqrt{1+\left(\frac{v_{s}}{v_m}\right)^{2}}-{\frac{v_s}{v_m}}\right]\right)^{-1},\\\\
&v_m=\frac{2b}{B_0}|\tau-\tau_f|，B_{0}=\frac{c_{d} K_{\mathrm{B}} T}{v_{s} b^{\alpha^{2}}}\\\\
\end{align*}
$$

 In the equation $\tau_f$ corresponds to the strong pinning force, while $\tau_c$ This corresponds to weak pinning. For different material systems, the composition of these two types of forces may vary.
For example, in some FCC crystal structures, the Peierls barrier can be considered a weak interaction, while dislocation junctions act as strong pinning obstacles.
The gradient form of the above equation is:

$$ \frac{\partial t_w}{\partial\tau}=-\frac{1}{\nu_0kT}\exp(\frac{Q_a}{kT})\frac{\partial Q_a}{\partial\tau}; $$

$$ \frac{\partial Q_a}{\partial\tau}=-\text{sgn}(\tau)\xi Q_0(\frac{||\tau|-\tau_c|}{\tau_o})^{\xi-1} $$

$$ \frac{\partial t_r}{\partial \tau}=-\text{sgn}(\tau)\frac{2bL}{v_m^2B_0}\frac{v_s^2}{v_m^2}(1-\frac{v_s}{v_m\sqrt{1+\left(\frac{v_s}{v_m}\right)^2}}) $$

### Some previously attempted but later discarded forms

The following forms consider only the weak or strong interactions individually, without establishing a consistent coupling between the weak and strong mechanisms.

$$
t_w=[\nu_D\frac {b^2 l}{\Omega_k}\exp(- Q_s/k_BT)]^{-1}\\ Q_s=(Q_0) [1-(|\tau| / \tau_c)^\alpha],Q_0>2k_BT\\  \Omega_k=\frac{Q_0}{\tau_P},\tau_c=\tau_P+\tau_b
$$

$$
t_{r}=\lambda \left(v_{s}\left[\sqrt{1+\left(\frac{v_{s}}{v_m}\right)^{2}}-{\frac{v_s}{v_m}}\right]\right)^{-1} \\
v_m=\frac{2b}{B_0}(|\tau|-\tau_c)+v_c，
B_{0}=\frac{c_{d} K_{\mathrm{B}} T}{v_{s} b^{\alpha^{2}}}
$$

$$
\frac{\partial t_w}{\partial\tau} = -\frac{sgn(\tau)\alpha Q_0^2}{\nu_D b^2 l \tau_P\tau_ck_BT}\exp(\frac{Q_s}{k_B T})(|\tau|/\tau_c)^{\alpha-1}
$$

$$
\frac{\partial t_r}{\partial\tau}=-sgn(\tau)\frac{2b\lambda}{v_m^2B_0}(\frac{v_s^2}{v_m^2}(1-\frac{v_s}{v_m}(1+(\frac{v_s}{v_m})^2)^{-1/2})
$$

or：

$$ t_w=[\nu_D\frac {b^2 l}{\Omega_k}\exp(- Q_s/k_BT)]^{-1} $$

$$ Q_s=(Q_0) [1-(\frac{|\tau|-\tau_b} {\tau_P})^\alpha],Q_0>2k_BT $$

$$ \Omega_k=\frac{Q_0}{\tau_P},\tau_c=\tau_P+\tau_b $$

The corresponding gradient changes to: 

$$
\frac{\partial t_w}{\partial\tau} = -\frac{\text{sgn}(\tau)\alpha Q_0^2}{\nu_D b^2 l \tau_P^2k_BT}\exp(\frac{Q_s}{k_B T})(\frac{|\tau|-\tau_b}{\tau_P})^{\alpha-1}
$$

# Hardening mechanism

$\rho_h$ represents the obstruction effect of dislocations from other slip systems on the current slip system, indicating mutual exclusivity； $\rho_J$  represents the contribution of the dislocation level on the current slip system to the interaction, allowing reduced obstructive effects from interactions when the dislocation density on the current system is low.

$$ \tau_c = \tau_P + c_b b G \sqrt{\rho_h + \rho_J} $$

$$ \rho_J = \sum_{\beta\neq\alpha} {h^{\alpha\beta} \sqrt{\rho^\beta_s\rho^\alpha_s}},\rho_s=\rho-\rho_0 $$

$$ \rho_h = \sum h^{\alpha\beta} \rho^\beta_e , \rho_e = \rho + \text{min} (\rho^\beta, \rho^\gamma), \alpha, \beta, \gamma \text{ are coplanar slips.} $$

To determine the type of each slip system, and to avoid initializing a massive matrix and manually mapping it to each slip system, it is preferable to identify the slip system behavior based on geometric criteria, $f^{\alpha\beta}$ classification and assignment of values：

1. Determine whether a dislocation reaction can occur spontaneously (Frank's rule) $b^\alpha\cdot b^\beta\ne0$ : if satisfied, continue; otherwise, apply H.
2. Determine whether the Burgers vectors are parallel $|\cos<b^\alpha,b^\beta>|=1$ ： If the vectors are parallel, assign to category N; otherwise, continue.
3. Determine whether two slip systems are coplanar by checking if $|\cos<n^\alpha,n^\beta>|=1$ ： If yes, assign category ： C， Otherwise, continue to the next check.
4. Determine whether the newly generated dislocation can slip on the original slip system $n^\alpha\cdot(b^\alpha+b^\beta)=0 || n^\beta\cdot(b^\alpha+b^\beta)=0$：If yes, assign G; otherwise, assign S.

Explanation of the above five interaction modes：

1. N, No junction: the Burgers vector of the new dislocation generated by the reaction of dislocations on two slip systems is parallel to the Burgers vector of the original slip system.
2. H, Hirth lock: the Burgers vector of the newly formed dislocation does not satisfy Frank’s rule (the new dislocation has higher energy than before the reaction).
3. C, Coplanar junction: the Burgers vector of the newly formed dislocation is coplanar with that of the original dislocation.
4. G, Glissile junction: the Burgers vector of the newly formed dislocation satisfies the energy criterion and can glide on the original slip plane.
5. S, Sessile junction: the Burgers vector of the newly formed dislocation satisfies the energy criterion, but it does not lie on either of the two slip planes of the original dislocations, and thus cannot glide.

For the above five cases, the cooperative hardening parameter correction factors are taken as: $a_N,a_{HL},a_{CJ},a_{GJ},a_{SJ}$ ，then we have： $a_{SJ}\gt a_{GJ} \gt a_{CJ} \gt a_{HL} \gt a_{N}$ 。

### Dependent on the saturated dislocation density

The saturated dislocation density in the material is rate-dependent, and according to[(Beyerlein & Tomé, 2008)](https://doi.org/10.1016/j.ijplas.2007.07.017) denoted as：

$$
\frac1{\sqrt{\rho_{sat}}}=\frac{k_2}{k_1}=\frac{c_h b}{g}(1-\frac{kT}{Db^3}\ln(\frac{|\dot\varepsilon|}{\dot\varepsilon_0}))
$$

However, this formulation still has details that need further discussion. The specific method for updating the dislocation density is also subject to various interpretations. In this code, the following form is adopted:

$$
d \rho=\left(k_{\rho, n u c} \frac{\left|\tau-\tau_{c, n u c}\right|}{G b^2}+\frac{k_{\rho, m u l}}{\bar{L}}\right)\left[1-\left(\frac{c_h b}{g}\left(1-\frac{k T}{D b^3} \log \frac{|\dot{\varepsilon}|}{\dot{\varepsilon}_0}\right)\right)^2 \rho\right] d \gamma
$$

The above equation includes a term for heterogeneous dislocation nucleation $k_{\rho, n u c} \frac{\left|\tau-\tau_{c, n u c}\right|}{G b^2}$ ，multiplication term $\frac{k_{\rho, m u l}}{\bar{L}}$ ，and control the final saturated dislocation density level.。
