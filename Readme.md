## Model Parameters and State Representation

Let:
- $N$ be the number of MS loci
- $r$ be the exponential tumor growth rate (per year)
- $\mu_i, \mu_d, \mu_s$ be the insertion, deletion, and SNV rates per MS locus per cell division, respectively
- $t_0, t_m, t_s$ be the initial time point, time when exponential growth begins, and patient's age at sampling, respectively
- $\Delta t$ be the time step for discrete simulations

The state of the tumor at time $t$ is represented by:
- $X(t)$: total number of cells
- $M_j(t) = [L_j(t), S_j(t)]$: state of MS locus $j$, where $L_j(t)$ is the length and $S_j(t)$ is the number of SNVs

## Model Dynamics

1. Population Growth:
   $$X(t + \Delta t) = X(t) \cdot e^{r \Delta t}$$

2. Mutation Events:
   For each MS locus $j$ (1 ≤ j ≤ N) in each time step $\Delta t$:
   $$\begin{aligned}
   p_i &= \mu_i \cdot \Delta t \\
   p_d &= \mu_d \cdot \Delta t \\
   p_s &= \mu_s \cdot \Delta t \\
   p_n &= 1 - (p_i + p_d + p_s)
   \end{aligned}$$
   $$(E_j) \sim \text{Multinomial}(X(t), [p_n, p_i, p_d, p_s])$$
   where $E_j = [E_{j,n}, E_{j,i}, E_{j,d}, E_{j,s}]$ represents the number of cells experiencing no event, insertion, deletion, and SNV respectively.

3. State Update:
   $$\begin{aligned}
   L_j(t + \Delta t) &= L_j(t) + E_{j,i} - E_{j,d} \\
   S_j(t + \Delta t) &= S_j(t) + E_{j,s}
   \end{aligned}$$

## State Distribution at $t = t_m$

To determine the state distribution at $t = t_m$, we define the following differential equations:
$$\begin{aligned}
\frac{dP_L}{dt} &= \mu_i P_{L-1} + \mu_d P_{L+1} - (\mu_i + \mu_d) P_L \\
\frac{dP_S}{dt} &= \mu_s P_{S-1} - \mu_s P_S
\end{aligned}$$

The probability distribution at $t_m$ is obtained by solving these equations:
$$P(M_j(t_m) = [L, S]) = f_j(L, S | \mu_i, \mu_d, \mu_s, t_m - t_0)$$

The initial state of each microsatellite locus at $t_m$ is drawn from this distribution:
$$M_j(t_m) \sim P(M_j(t_m) = [L, S])$$

## Stochastic Simulation Algorithm

1. Initialize $X(t_0)$ and $M_j(t_0)$ for all $j$
2. For $t = t_0$ to $t_s$ in steps of $\Delta t$:
   1. If $t \geq t_m$: Update $X(t + \Delta t) = X(t) \cdot e^{r \Delta t}$
   2. For each MS locus $j$:
      - Generate $E_j$ from Multinomial distribution
      - Update $M_j(t + \Delta t)$
   3. Record $X(t)$ and $M_j(t)$ for all $j$

## Probability Distribution Estimation

To estimate the probability distribution of MS loci states:
1. Generate $K$ independent stochastic paths using the above algorithm.
2. For each locus $j$, create a 2D histogram of $[L_j, S_j]$ values across all $K$ paths.
3. Normalize the histogram to obtain an estimated probability mass function (PMF) $P_j(L, S)$.

## Likelihood Function

Let $O_j$ be the observed state of locus $j$ in a tumor sample. The likelihood function is:
$$L(\theta | O) = \prod_{j=1}^N P_j(O_j | \theta)$$
where $\theta$ represents the model parameters $\{r, \mu_i, \mu_d, \mu_s, t_m\}$.

The log-likelihood is:
$$\log L(\theta | O) = \sum_{j=1}^N \log P_j(O_j | \theta)$$

## Bias Correction

To account for the stochastic nature of the simulation, we apply a bias correction:
$$\log L_{\text{corrected}}(\theta | O) = \log L(\theta | O) + \frac{c}{K}$$
where $c$ is a constant estimated empirically.

## Sample Purity Consideration

To account for normal cell contamination, we introduce a purity parameter $\rho$ and use a probabilistic approach:
$$P(O_j^{\text{measured}} | \rho, \theta) = \rho P(O_j^{\text{tumor}} | \theta) + (1-\rho) P(O_j^{\text{normal}})$$
where $P(O_j^{\text{tumor}} | \theta)$ is the probability of observing the tumor state given the model parameters, and $P(O_j^{\text{normal}})$ is the probability of the normal state.

## Extended Model for Testing MS Equilibrium

### Equilibrium Dynamics

Under equilibrium, the mutation rates satisfy:
$$\mu_i P(L - 1) = \mu_d P(L + 1)$$
for insertions and deletions at length $L$, and:
$$\mu_s P(S - 1) = \mu_s P(S + 1)$$
for SNVs.

### Likelihood Function for Model Comparison

1. Equilibrium Model:
   $$L_{\text{eq}}(\theta | O) = \prod_{j=1}^N P_{\text{eq}}(O_j | \theta)$$

2. Non-Equilibrium Model:
   $$L_{\text{non-eq}}(\theta | O) = \prod_{j=1}^N P_{\text{non-eq}}(O_j | \theta)$$

3. Log-Likelihood Comparison:
   - Equilibrium: $LL_{\text{eq}} = \log L_{\text{eq}}(\theta | O)$
   - Non-Equilibrium: $LL_{\text{non-eq}} = \log L_{\text{non-eq}}(\theta | O) + c/K$

### Statistical Test for Equilibrium

Likelihood Ratio Test (LRT):
$$LR = 2 (LL_{\text{non-eq}} - LL_{\text{eq}})$$

Under the null hypothesis (equilibrium), $LR \sim \chi^2_k$, where $k = d_f(\text{non-eq}) - d_f(\text{eq})$.
