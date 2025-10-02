
# Quantile Regression for Cancer Research: Beyond the Mean


> **Quantile Regression in Cancer Research: Modeling Effects Across the Outcome Distribution**

Quantile regression helps you understand how covariates relate to **different parts of an outcome distribution** (e.g., the **lower tail** vs the **upper tail**) rather than just the mean. This is crucial in cancer studies where clinical interest often centers on **high-risk** or **vulnerable** subgroups (e.g., patients with **short telomeres**, **very low BMI**, or **long survival**).



## Why not just model the mean?

Classical regression targets the conditional mean:  
$E(Y\mid X)=X\beta$

This answers: “On average, how does $Y$ change with $X$?” But it **does not** tell you how $X$ relates to:
- Patients with **unusually low** outcomes (e.g., **lower BMI**, **shorter telomere length**),
- Patients with **unusually high** outcomes (e.g., **higher BMI**, **longer telomere length**),
- Or how effects **vary across the distribution** (heterogeneous effects).

When effects are **not constant** across the distribution (common in biomedical data), mean regression can **mask clinically important patterns**.



## What is quantile regression?

Quantile regression models a **conditional quantile** of $Y$ given $X$:

$$
Q_Y(\tau\mid X)=X\,\beta(\tau),\quad \tau\in(0,1)
$$

- $Q_Y(0.5\mid X)$ is the **conditional median**.  
- $Q_Y(0.1\mid X)$ targets the **lower 10th percentile** (e.g., “lower tail”).  
- $Q_Y(0.9\mid X)$ targets the **upper 90th percentile** (e.g., “upper tail”).

Because $\beta(\tau)$ **depends on $\tau$**, the effect of a covariate can differ at the lower vs upper tails—exactly the insight mean models miss.

**Interpretation (mean vs quantile models)**

- **Mean model (ordinary regression)**
  - Targets the **average** outcome: $E(Y\mid X)=X\beta$.
  - A coefficient $\beta_j=0.5$ means: **for a 1-unit increase in $X_j$, the mean of $Y$ increases by 0.5 units**, holding other variables fixed.
  - Example (BMI): If `male` has $\beta_{\text{male}}=+0.4$, then **on average** males have BMI higher by 0.4 units than females with the same covariates.

- **Quantile model (regression at a specific percentile)**
  - Targets a **percentile** of the outcome: $Q_Y(\tau\mid X)=X\beta(\tau)$.
  - A coefficient $\beta_j(\tau)$ tells how a 1-unit change in $X_j$ **moves the $\tau$-th percentile** of $Y$ (not the mean), with other variables fixed.

  - **Lower tail (e.g., $\tau=0.10$):**  
    If $\beta_j(0.10)=+0.8$, then a 1-unit increase in $X_j$ is associated with the **10th percentile** of $Y$ being **0.8 units higher** (a shift **up** in the lower end).

  - **Upper tail (e.g., $\tau=0.90$):**  
    If $\beta_j(0.90)=-0.2$, then a 1-unit increase in $X_j$ is associated with the **90th percentile** of $Y$ being **0.2 units lower** (a shift **down** in the upper end).

- **Heterogeneous effects across the distribution**
  - If $\beta_j(0.10)$ and $\beta_j(0.90)$ differ in **size or sign**, the effect of $X_j$ is **distributionally heterogeneous**—it impacts low and high values of $Y$ differently.
  - Example (telomere length): `smoking` shows $\beta_{\text{smk}}(0.10)=-1.5$ but $\beta_{\text{smk}}(0.90)=-0.2$. Interpretation: smoking is **much more detrimental among the shortest telomere lengths** than near the longest.




## More Examples in Cancer Research

1. **Telomere length (continuous)**  
   *Is the effect of sex different at the **shortest** versus **median** telomere lengths?*  
   Quantile regression can reveal whether sex differences are larger among **short-telomere** patients (e.g., $\tau=0.10$) than around the median ($\tau=0.50$).

2. **Identifying CpG Sites Associated with High LINE-1 Activity (via total TE counts)**

**Context.** LINE-1 activity is **not directly observed**. We use **total transposable-element (TE) counts** as a *proxy* outcome to infer patterns consistent with **high LINE-1 activity**, recognizing that total TE counts aggregate LINE-1 and non-LINE-1 families.

**Working assumption.** Increases in LINE-1 activity tend to **raise** total TE counts (on average), so associations that appear **strongest in the upper tail** of total TE counts are informative about periods of **elevated LINE-1 activity**.

**Question.** Are methylation levels at specific CpG sites associated with the **upper quantiles** of **total TE counts** (e.g., $\tau=0.80$–$0.95$)?



3. **Survival times (with censoring)**  
   *Do race/ethnicity effects differ for **short** versus **long** survival?*  
   Standard quantile regression assumes fully observed outcomes; for right-censored survival, use **censored quantile regression** or other survival-adapted quantile methods. The idea—**effects vary across survival quantiles**—still applies, but the model must handle censoring.



## R Example

```r
# install.packages("quantreg")  # if needed
library(quantreg)

# Example: BMI ~ sex + age + stage + smoking
taus <- c(0.10, 0.50, 0.90)

fits <- lapply(taus, function(tau) {
  rq(BMI ~ sex + age + stage + smoking, tau = tau, data = dat)
})

# Summaries
lapply(fits, summary)

# Compare the sex coefficient across quantiles
sapply(fits, function(f) coef(f)["sexMale"])
````

> **Tip:** Visualize coefficients over many $\tau$ values (e.g., 0.05–0.95) to see where effects are strongest. In R, loop over `seq(0.05, 0.95, by = 0.05)` and plot the coefficient paths.



## References (starter list)

* Koenker, R., & Bassett, G. (1978). **Regression quantiles.** *Econometrica*, 46(1), 33–50.
* Koenker, R. (2005). **Quantile Regression.** Cambridge University Press.
* Hong, H. G., Christiani, D. C., & Li, Y. (2019). *Quantile regression for survival data in modern cancer research: expanding statistical tools for precision medicine.* Precision Clinical Medicine, 2(3), 90–99. [PMCID: PMC6644129](https://pmc.ncbi.nlm.nih.gov/articles/PMC6644129/)

