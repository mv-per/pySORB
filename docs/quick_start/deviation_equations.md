# Deviation Equations

## Supported Deviation Equations

The following deviation equations are currently supported[^1]:

- Sum of Squared Errors (`'SSE'`)

$$\sigma = \sum_{i}^{N} \dfrac{Calc[i] - Exp[i]}{2.0}^2$$

- Average Relative Error (`'ARE'`)

$$\sigma = \frac{100.0}{N} \times \sum_{i}^{N} \left(\frac{Exp[i] - Calc[i]}{Exp[i]}\right)$$

- Absolute Error (`'EABS'`)

$$\sigma = \sum_{i}^{N} |Exp[i] - Calc[i]|$$

- Relative Absolute Error (`'RABS'`)

$$\sigma = \sum_{i}^{N} \times \dfrac{\left| Exp[i] - Calc[i] \right|}{Exp[i]}$$

- Hybrid Deviation (`'HYBRID'`)

$$\sigma = \dfrac{100}{N-K} \times \sum_{i}^{N} \dfrac{Exp[i] - Calc[i]}{Exp[i]}^2$$

- Mean Percentage Standard Deviation (`'MPSD'`)

$$\sigma = 100 \sqrt{\sum_{i=0}^{N-1} \left(\frac{{(Exp[i] - Calc[i])^2}}{{Exp[i]}}\right)}$$

- Sum of Relative Errors (`'SRE'`)

$$\sigma_{ARE} = \frac{100.0}{N} \times \sum_{i}^{N} \left(\frac{Exp[i] - Calc[i]}{Exp[i]}\right)$$
$$\sigma = \dfrac{\sqrt{\sum_{i}^{N} \left(\frac{{(Exp[i] - Calc[i])}}{{Exp[i]}} \sigma_{ARE}\right)^2}}{\sqrt{N-1}}$$

- Chi-square Deviation ('`CHI_2'`)

$$\sigma = \sum_{i}^{N} \frac{{(Exp[i] - Calc[i])^2}}{{Exp[i]}}$$

- Residual Sum of Squares Deviation (`'R_S'`)

$$\sigma_{sse} = \sum_{i}^{N} \dfrac{Calc[i] - Exp[i]}{2.0}^2$$
$$\sigma = 1.0 - \frac{{6.0 \sigma_{sse}}}{N} \cdot \left(N - 1.0\right)^2$$

## References

[^1]: Foo, K. Y., & Hameed, B. H. (2010). Insights into the modeling of adsorption isotherm systems. Chemical Engineering Journal, 156(1), 2–10. doi:10.1016/j.cej.2009.09.013
