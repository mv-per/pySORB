# Fitting

n (mol/kg) = f(P[Pa], T[K])

## Empirical isotherms

The following classic adsorption isotherm models are currently supported:

|                     Isotherm                     |          Code          |                      Parameters                       | Equation                                                                                                                                                                                       |
| :----------------------------------------------: | :--------------------: | :---------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|                    Freundlich                    |     `'freundlich'`     |                       $K$, $t$                        | $n = KP^{1/t}$                                                                                                                                                                                 |
| Temperature-dependent Freundlich[^DUONG-DO-BOOK] |    `'freundlich-2'`    |                $K_0$, $\alpha$, $A_0$                 | $a = \dfrac{RT}{A_0}$ <br> $K = K_0 e^{-\alpha f_n}$ <br> $ n = KP^{a}$                                                                                                                        |
|                     Langmuir                     |      `'langmuir'`      |                    $n_{max}$, $b$                     | $n = n_{max} \dfrac{bP}{(1+bP)}$                                                                                                                                                               |
|  Temperature-dependent Langmuir[^DUONG-DO-BOOK]  |     `'langmuir-2'`     |             $n_{max}$, $b_{\infty}$, $Q$              | $b = b_\infty e^{Q/RT}$ <br> $n = n_{max} \dfrac{bP}{(1+bP)}$                                                                                                                                  |
|                  Dual-Langmuir                   |   `'dual-langmuir'`    |        $n_{max,1}$, $b_1$ , $n_{max,2}$, $b_2$        | $n = n_{max,1}  \dfrac{b_1P}{(1+b_1P)} + n_{max,2}\dfrac{b_2P}{(1+b_2P)}$                                                                                                                      |
|                       Sips                       |        `'sips'`        |                 $n_{max}$, $b$ , $t$                  | $n = n\_{max} \dfrac{(bP)^{1/t}}{1+(bP)^{1/t}} $                                                                                                                                               |
|    Temperature-depentent Sips[^DUONG-DO-BOOK]    |       `'sips-2'`       | $n_{max,\infty}$, $\chi$, $t_0$, $\alpha$, $b_0$, $Q$ | $T_0 = 298.15$ <br> $n_{max} = n_{max,\infty}e^{[\chi(1-T/T_0)]}$ <br> $b = b_0 e^{[Q/RT(T_0/T-1)]}$ <br> $1/t = 1/t_0 + \alpha(1-T_0/T)$ <br>$n = n\_{max} \dfrac{(bP)^{1/t}}{1+(bP)^{1/t}} $ |
|                       Toth                       |        `'toth'`        |                 $n_{max}$, $b$ , $t$                  | $n = n\_{max} \dfrac{bP}{(1+(bP)^t)^\tfrac{1}{t}} $                                                                                                                                            |
|    Temperature-dependent Toth[^DUONG-DO-BOOK]    |       `'toth-2'`       | $n_{max,\infty}$, $\chi$, $t_0$, $\alpha$, $b_0$, $Q$ | $T_0 = 298.15$ <br> $n_{max} = n_{max,\infty}e^{[\chi(1-T/T_0)]}$ <br> $b = b_0 e^{[Q/RT(T_0/T-1)]}$ <br> $t = t_0 + \alpha(1-T_0/T)$ <br> $n = n\_{max} \dfrac{bP}{(1+(bP)^t)^\tfrac{1}{t}} $ |
|                      Unilan                      |       `'unilan'`       |               $n_{max}$, $\bar{b}$, $t$               | $n = \dfrac{n_{max}}{2t} ln\left(\dfrac{1+\bar{b}e^tP}{1+\bar{b}e^{-t}P} \right)$                                                                                                              |
|                Keller-Staudt-Toth                | `'keller-staudt-toth'` |         $n_{max}$, $b$, $\alpha_m$, $\beta_0$         | $\alpha = \dfrac{1+\alpha_m\beta_0}{1+\beta_0 P}$<br> $n = n_{max}\dfrac{\alpha_mbP}{(1+(bP)^\alpha)^{1/\alpha}}$                                                                              |
|                 Redlich-Peterson                 | `'redlich-perterson'`  |                     $K$, $b$, $t$                     | $n = \dfrac{KP}{(1+(bP)^t)} $                                                                                                                                                                  |
|                  Jensen-Seaton                   |   `'jensen-seaton'`    |                $K$, $a$, $\kappa$, $t$                | $n = KP \left(1 + \left[\dfrac{KP}{a(1+\kappa P)}\right]^t  \right)^{-1/t}$                                                                                                                    |

Parameter Units

|              Parameter              |         Unit          |
| :---------------------------------: | :-------------------: |
|           $n_{max}$, $a$            |     $mol.kg^{-1}$     |
| $b$, $\bar{b}$, $\kappa$, $\beta_0$ |       $Pa^{-1}$       |
|                 $K$                 | $mol.kg^{-1}.Pa^{-1}$ |
|           $t$, $\alpha_m$           |           -           |
|                  Q                  |     $J.mol^{-1}$      |

## Vacancy Solution Method

$$
\theta = \dfrac{n}{n_{max}}
$$

Activity coefficients:

- Wilson (`'wilson'`)[^VSM-WILSON-PURE]:

  Optimization parameters: $n_{max}$, $b$, $\alpha_{v1}$, $\alpha_{1v}$.

  $$
  p = \left[ \dfrac{n_{max}}{b} \dfrac{\theta}{1-\theta} \right] \left[ \alpha_{1v} \dfrac{1-(1-\alpha_{v1})\theta}{\alpha_{1v} + (1-\alpha_{1v})\theta} \right] \exp\left[-\dfrac{\alpha_{v1}(1-\alpha_{v1})\theta}{1-(1-\alpha{v1})\theta} - \dfrac{(1-\alpha_{1v})\theta}{\alpha_{1v}+(1-\alpha_{1v})\theta} \right]
  $$

- NRTL (`'nrtl'`)[^VSM-NRTL-PURE]:

  Optimization parameters: $n_{max}$, $b$, $\alpha_{1v}$, $\tau_{1v}$, $\tau_{v1}$.

  $$
  G_{1v} = \exp(-\alpha_{1v}\tau_{1v})\\
  G_{v1} = \exp(-\alpha_{1v}\tau_{v1})\\
  $$

  $$
    p = \left[ \dfrac{n_{max}}{b} \dfrac{\theta}{1-\theta} \right] \exp \left[ \dfrac{\tau_{1v}G_{1v}^2}{(G_{1v}-1){(G_{1v}-1)\theta+1}^2} - \dfrac{\tau_{v1}G_{v1}^2}{(G_{1v}-1){(G_{v1}-1)\theta-G\_{v1}}^2} - \dfrac{\tau_{v1}(G_{1v}-1)+\tau_{1v}G_{1v}^2(G_{v1}-1)}{(G*{1v}-1)(G\_{v1}-1)}\right]
  $$

- Flory-Huggins (`'flory-huggins'`)[^VSM-FH-PURE]:

  Optimization parameters: $n_{max}$, $b$, $\alpha_{1v}$.

  $$
  p = \left[ \dfrac{n_{max}}{b} \dfrac{\theta}{1-\theta} \right] \exp\left[\dfrac{\alpha_{1v}^2\theta}{1+\alpha_{1v}\theta} \right]
  $$

## Potential Theory Models[^MPTA-ORIGINAL]

The Multicomponent potential theory of adsorption is based on the equilibrium of chemical potentials between the bulk and adsorbed phases, taking into account an adsorption potential that varies within the adsorbed phase.

$$
\mu_{ad} - \varepsilon_z = \mu_{b}
$$

Applying the iso-fugacity rule, we have:

$$
f_{ad}(P_z) = f_{b}(P_b)exp(f(\varepsilon_z))
$$

The above equation is solved in all the adsorbed phase volume (z) and the adsorbed loading is calculated integrating the density of adsorbate:

$$
n^{ex} = \int_0^z \rho_{ad}(f_{ad}(P_z)) - \rho_{b}(f_{b}(P_b))
$$

$$
n^{abs} = \int_0^z \rho_{ad}(f_{ad}(P_z))
$$

So, for each step in the integral, we must calculate the adsorption potential.

The pySORB currently supports the following adsorption potentials:

- Dubinin-Radushkevich-Astakhov (`'DRA`')[^DRA-POTENTIAL]

Optimization parameters: $\varepsilon_0$, $z_0$, $\beta$.

$$
\varepsilon = \varepsilon_0 ln \ \left( \dfrac{z_0}{z}\right)^{\dfrac{1}{\beta}}
$$

$$
f(\varepsilon_z) = \dfrac{\varepsilon}{RT}
$$

- LEE (`'LEE`')[^LEE-POTENTIAL]:

Parameters = $\varepsilon_{fs}$, $L$, $A$

The Lee potential is based on the Lennard-Jonnes potential between two solid walls and has a shape similar to:

$\sigma_{fs}$: fluid-solid colision diameter (Angstrom)

$$
\sigma_{fs} = \dfrac{\sigma_{ff} + \sigma_{ss}}{2}
$$

$\sigma_{ff}$: Adsorbent wall lennard-jonnes diameter
$\sigma_{ss}$: Adsorbate lennard-jonnes Diameter

$$
\epsilon_{fs} = \varepsilon_{fs}/k_b
$$

$k_b$: Boltzmann's constant

$$
\psi = 4 \pi \rho_{atoms} \epsilon_{fs} \sigma_{fs}^2 \left(\dfrac{\sigma_{fs}^{10}}{5z^{'10}}-\dfrac{1}{2}\sum_{i=0}^{3}\dfrac{\sigma_{fs}^4}{[z^{'} + i\sigma_{ss}]^4} \right)
$$

$$
f(\varepsilon_z) = \dfrac{\psi_z + \psi_{L-z}}{k_{B}T}
$$

$$
n^{ex} = A 2.0 \int_{0.7\sigma_{ss}}^\dfrac{L}{2} \rho_{ad}(f_{ad}(P_z)) - \rho_{b}(f_{b}(P_b))
$$

$$
n^{abs} = A 2.0 \int_{0.7\sigma_{ss}}^\dfrac{L}{2} \rho_{ad}(f_{ad}(P_z))
$$

- STEELE-10-4-3 (`'STEELE`')[^STEELE-POTENTIAL]:

Parameters = $\varepsilon_{fs}$, $L$, $A$

$$
\epsilon_{fs} = \varepsilon_{fs}/k_b
$$

$k_b$: Boltzmann's constant

$$
\psi = 2\pi\rho_s\epsilon_{fs}\sigma_{fs}^2\Delta \left(\dfrac{2}{5} \left[ \dfrac{\sigma_{fs}}{z}\right]^{10} - \left[\dfrac{\sigma_{fs}}{z} \right]^{4} - \left[\dfrac{\sigma_{fs}^4}{3-\Delta(z+0.61\Delta)^3} \right] \right)
$$

$$
f(\varepsilon_z) = \dfrac{\psi_z + \psi_{L-z}}{k_{B}T}
$$

$$
n^{ex} = A 2.0 \int_{0.7\sigma_{ss}}^\dfrac{L}{2} \rho_{ad}(f_{ad}(P_z)) - \rho_{b}(f_{b}(P_b))
$$

$$
n^{abs} = A 2.0 \int_{0.7\sigma_{ss}}^\dfrac{L}{2} \rho_{ad}(f_{ad}(P_z))
$$

[^DUONG-DO-BOOK]: Do, Duong D. Adsorption analysis: Equilibria and kinetics. Vol. 2. World Scientific, 1998.
[^MPTA-ORIGINAL]: Shapiro, A. A., & Stenby, E. H. (1998). Potential Theory of Multicomponent Adsorption. Journal of Colloid and Interface Science, 201(2), 146–157. doi:10.1006/jcis.1998.5424
[^VSM-WILSON-PURE]: Suwanayuen, S., & Danner, R. P. (1980). A gas adsorption isotherm equation based on vacancy solution theory. AIChE Journal, 26(1), 68–76. doi:10.1002/aic.690260112
[^VSM-WILSON-MIX]: Suwanayuen, S., & Danner, R. P. (1980). Vacancy solution theory of adsorption from gas mixtures. AIChE Journal, 26(1), 76–83. doi:10.1002/aic.690260113
[^VSM-NRTL-PURE]: Munakata, K. (2007). Vacancy Solution Model Formulated by the NRTL Equation for Correlation of Adsorption Equilibria. Journal of Chemical Engineering of Japan, 40(5), 398–409. doi:10.1252/jcej.40.398
[^VSM-FH-PURE]: Cochran, T. W., Kabel, R. L., & Danner, R. P. (1985). Vacancy solution theory of adsorption using Flory-Huggins activity coefficient equations. AIChE Journal, 31(2), 268–277. doi:10.1002/aic.690310214
[^STEELE-POTENTIAL]: Steele, W. A. (1973). The physical interaction of gases with crystalline solids. Surface Science, 39(1), 149–175. doi:10.1016/0039-6028(73)90102-7
[^LEE-POTENTIAL]: Lee, Lloyd L. Molecular thermodynamics of nonideal fluids. Butterworth-Heinemann, 1988.
[^DRA-POTENTIAL]: Dubinin, M. M. (1983). Microporous structures and absorption properties of carbonaceous adsorbents. Carbon, 21(4), 359–366. doi:10.1016/0008-6223(83)90128-8.
