# Comprehensive Mathematical Model for Gold Scalping Strategies

> **Pure mathematical formulas, equations, and derivations for intraday gold scalping.**
> All notation follows standard financial mathematics and LaTeX conventions.

---

## Table of Contents

1. [Microstructure Mathematics](#1-microstructure-mathematics)
2. [High-Frequency Price Dynamics](#2-high-frequency-price-dynamics)
3. [Entry/Exit Optimization](#3-entryexit-optimization)
4. [Execution Mathematics](#4-execution-mathematics)
5. [Profitability Analysis](#5-profitability-analysis)
6. [Risk Management for Scalping](#6-risk-management-for-scalping)
7. [Volatility Models for Scalping](#7-volatility-models-for-scalping)

---

## 1. Microstructure Mathematics

### 1.1 Bid-Ask Spread Analysis

The **quoted spread** at time $t$ is defined as:

$$S(t) = A(t) - B(t)$$

where $A(t)$ is the best ask price and $B(t)$ is the best bid price.

The **relative (percentage) spread** normalised by the midpoint $M(t) = \tfrac{A(t)+B(t)}{2}$:

$$s(t) = \frac{S(t)}{M(t)} = \frac{A(t) - B(t)}{\tfrac{1}{2}[A(t) + B(t)]}$$

The **Roll (1984) implied spread** estimated from autocovariance of price changes $\Delta P_t$:

$$\hat{S} = 2\sqrt{-\text{Cov}(\Delta P_t,\, \Delta P_{t-1})}$$

**Spread decomposition** into adverse selection ($\lambda$), inventory ($\phi$), and order processing ($c$) components:

$$S(t) = \underbrace{2\lambda}_{\text{adverse selection}} + \underbrace{2\phi}_{\text{inventory}} + \underbrace{2c}_{\text{processing}}$$

**Expected round-trip cost** for a scalp trade of size $Q$:

$$C_{\text{spread}}(Q) = \frac{1}{2} S(t) \cdot Q$$

---

### 1.2 Order Flow Dynamics

Define the **signed order-flow** at tick $k$ as:

$$x_k = d_k \cdot V_k$$

where $d_k \in \{+1, -1\}$ is the trade direction (buyer/seller initiated) and $V_k$ is the traded volume.

**Cumulative order-flow** over the interval $[0,T]$:

$$\mathcal{X}(T) = \sum_{k=1}^{N(T)} x_k$$

**Order-flow imbalance** (OFI), a predictor of short-horizon price movement:

$$\text{OFI}(t) = \frac{\sum_{k} d_k V_k}{\sum_{k} V_k}$$

**Order-flow autocorrelation** at lag $h$:

$$\rho(h) = \frac{\mathbb{E}[x_k \, x_{k-h}] - \mu_x^2}{\sigma_x^2}$$

**Price impact of order flow** (Kyle 1985 linear model):

$$\Delta M(t) = \lambda \, \mathcal{X}(t) + \epsilon(t), \quad \epsilon(t) \sim \mathcal{N}(0, \sigma_\epsilon^2)$$

where $\lambda$ is Kyle's **lambda** (price impact coefficient):

$$\lambda = \frac{\text{Cov}(\Delta M, \mathcal{X})}{\text{Var}(\mathcal{X})}$$

---

### 1.3 Liquidity Models

**Amihud (2002) illiquidity ratio** for day $d$:

$$\text{ILLIQ}_d = \frac{|\Delta P_d|}{V_d}$$

**Average illiquidity** over $N$ trading days:

$$\overline{\text{ILLIQ}} = \frac{1}{N} \sum_{d=1}^{N} \text{ILLIQ}_d$$

**Liquidity-adjusted return**:

$$\tilde{r}_d = r_d - \frac{1}{2} \Delta \text{ILLIQ}_d \cdot W$$

where $W$ is the position size in currency units.

**Resilience function** — the rate at which the spread returns to equilibrium after a shock:

$$\frac{dS(t)}{dt} = -\kappa \bigl[S(t) - \bar{S}\bigr] + \eta(t), \quad \eta(t) \sim \mathcal{N}(0, \sigma_\eta^2)$$

with mean-reversion speed $\kappa > 0$ and long-run spread $\bar{S}$.

---

### 1.4 Market Depth Functions

**Aggregate depth** at price level $p$ on the bid and ask sides:

$$D^B(p, t) = \sum_{i} q_i^B \cdot \mathbf{1}\{P_i^B \geq p\}, \qquad D^A(p, t) = \sum_{j} q_j^A \cdot \mathbf{1}\{P_j^A \leq p\}$$

**Depth imbalance** at the best quotes:

$$\text{DI}(t) = \frac{D^B(B(t), t) - D^A(A(t), t)}{D^B(B(t), t) + D^A(A(t), t)}$$

**Limit order book (LOB) density** — exponential decay model for depth at distance $\delta$ from mid:

$$\rho(\delta) = \rho_0 \, e^{-\alpha \delta}, \qquad \alpha > 0$$

**Expected market-impact cost** of a market order of size $Q$ against the LOB:

$$C_{\text{MI}}(Q) = \int_0^Q \bigl[p(\xi) - M(t)\bigr] \, d\xi$$

where $p(\xi)$ is the marginal price obtained by consuming $\xi$ units of depth.

---

## 2. High-Frequency Price Dynamics

### 2.1 Tick-by-Tick Price Movement

**Continuous-time Brownian motion model** for the log-price $\ln P(t)$:

$$dP(t) = \mu \, P(t) \, dt + \sigma \, P(t) \, dW(t)$$

or equivalently for the log-price $X(t) = \ln P(t)$:

$$dX(t) = \left(\mu - \frac{\sigma^2}{2}\right) dt + \sigma \, dW(t)$$

where $W(t)$ is a standard Brownian motion, $\mu$ is drift, and $\sigma$ is diffusion coefficient.

**Discrete tick-by-tick representation** over interval $\Delta t$:

$$P_{n+1} = P_n + \mu \Delta t + \sigma \sqrt{\Delta t}\, Z_n, \quad Z_n \overset{\text{i.i.d.}}{\sim} \mathcal{N}(0,1)$$

**Realised variance** estimated from $N$ tick returns $r_i = \ln(P_{i}/P_{i-1})$:

$$\widehat{\sigma^2} = \text{RV} = \sum_{i=1}^{N} r_i^2$$

---

### 2.2 Jump Processes

**Jump-diffusion model** (Merton 1976):

$$dP(t) = \mu P(t)\,dt + \sigma P(t)\,dW(t) + P(t^-)\,dJ(t)$$

where $J(t) = \sum_{k=1}^{N(t)} (e^{Y_k} - 1)$ with:
- $N(t)$: Poisson process with intensity $\lambda$ (jumps per unit time)
- $Y_k \overset{\text{i.i.d.}}{\sim} \mathcal{N}(\mu_J, \sigma_J^2)$: i.i.d. log-jump sizes

**Compensated jump-diffusion** (zero mean):

$$dX(t) = \left(\mu - \frac{\sigma^2}{2} - \lambda \bar{Y}\right)dt + \sigma\,dW(t) + \sum_{k=1}^{dN(t)} Y_k$$

**Jump intensity** estimated from intraday data (Lee–Mykland 2008 test statistic):

$$\mathcal{L}(i) = \frac{|r_i|}{\hat{\sigma}(i)}, \qquad \hat{\sigma}(i) = \frac{1}{\sqrt{K-2}} \sum_{j=i-K+1}^{i-2} |r_j|$$

A jump is detected at tick $i$ if $\mathcal{L}(i) > c_\alpha$ where $c_\alpha$ is the critical value at significance level $\alpha$.

**Expected number of jumps** in time horizon $T$:

$$\mathbb{E}[N(T)] = \lambda T$$

---

### 2.3 Ultra-High Frequency Models

**Autoregressive Conditional Duration (ACD) model** — models the time $\tau_n$ between trades:

$$\tau_n = \psi_n \, \varepsilon_n, \qquad \varepsilon_n \overset{\text{i.i.d.}}{\sim} \text{Exp}(1)$$

$$\psi_n = \omega + \alpha \tau_{n-1} + \beta \psi_{n-1}$$

with $\omega > 0$, $\alpha \geq 0$, $\beta \geq 0$, and $\alpha + \beta < 1$ for stationarity.

**Hawkes process** — self-exciting point process for trade arrivals with intensity:

$$\lambda(t) = \mu + \int_{-\infty}^{t} \phi(t - s) \, dN(s)$$

with exponential kernel $\phi(u) = \alpha e^{-\beta u}$, $\alpha < \beta$.

**Long-run (unconditional) intensity**:

$$\bar{\lambda} = \frac{\mu}{1 - \alpha/\beta}$$

---

### 2.4 Price Impact Functions

**Power-law temporary price impact** of a market order of signed volume $Q$:

$$I(Q) = \alpha \, |Q|^\beta \cdot \text{sgn}(Q), \qquad 0 < \beta \leq 1$$

For gold futures, empirical estimates: $\beta \approx 0.5$ (square-root impact law).

**Almgren–Chriss linear temporary impact**:

$$I_{\text{temp}}(v) = \eta \, v + \gamma \frac{v^2}{2}$$

where $v = Q/T$ is the trading rate, $\eta$ is temporary impact coefficient, $\gamma$ is permanent impact coefficient.

**Total execution cost** for a trajectory liquidating $X_0$ shares over $[0,T]$:

$$C(v) = \int_0^T \left[\eta\, v(t)^2 + \gamma \, v(t) \, X(t)\right] dt$$

**Square-root (Gabaix et al.) impact**:

$$\frac{\Delta P}{P} = Y \cdot \sqrt{\frac{Q}{V_{\text{day}}}}$$

where $Y \approx 0.1$–$0.3$ for most markets and $V_{\text{day}}$ is the average daily volume.

---

## 3. Entry/Exit Optimization

### 3.1 Optimal Entry Points using Bellman Equations

Let $V(x, t)$ be the value function of the scalper's optimisation problem. The **Hamilton–Jacobi–Bellman (HJB) equation** is:

$$\frac{\partial V}{\partial t} + \sup_{u \in \mathcal{U}} \left\{ \mathcal{L}^u V(x,t) + r(x, u) \right\} = 0$$

with terminal condition $V(x, T) = g(x)$, where:
- $\mathcal{L}^u$ is the infinitesimal generator under control $u$
- $r(x, u)$ is the instantaneous reward
- $g(x)$ is the terminal payoff

**Verification theorem**: If $V \in C^{1,2}$ solves the HJB equation, then the optimal control is:

$$u^*(x,t) = \arg\sup_{u \in \mathcal{U}} \left\{ \mathcal{L}^u V(x,t) + r(x,u) \right\}$$

**Optimal entry threshold** $P^*_{\text{entry}}$ satisfies:

$$V(P^*_{\text{entry}},\, t) = h(P^*_{\text{entry}})$$

$$V'(P^*_{\text{entry}},\, t) = h'(P^*_{\text{entry}}) \quad \text{(smooth-pasting condition)}$$

where $h(P)$ is the value of entering the trade immediately.

---

### 3.2 Optimal Exit Rules

**Profit-taking threshold** $P^*_{\text{TP}}$ and **stop-loss threshold** $P^*_{\text{SL}}$ solve:

$$P^*_{\text{TP}} = \arg\max_{P > P_{\text{entry}}} \mathbb{E}\left[e^{-r\tau}(P_\tau - P_{\text{entry}} - C)\right]$$

$$P^*_{\text{SL}} = \arg\max_{P < P_{\text{entry}}} \mathbb{E}\left[e^{-r\tau}(P_{\text{entry}} - P_\tau - C)\right]$$

where $C$ is the round-trip transaction cost and $\tau$ is the stopping time.

**Optimal stopping time**:

$$\tau^* = \inf \left\{ t \geq t_0 : P(t) \geq P^*_{\text{TP}} \text{ or } P(t) \leq P^*_{\text{SL}} \right\}$$

**Expected trade duration** (GBM with two-sided exit):

Let $p_{\text{TP}} = P(P_\tau = P^*_{\text{TP}})$ be the probability of exiting at the take-profit barrier and $p_{\text{SL}} = 1 - p_{\text{TP}}$ be the probability of exiting at the stop-loss barrier. For a GBM with drift $\mu$ and volatility $\sigma$:

$$p_{\text{TP}} = \frac{(P_0 / P^*_{\text{SL}})^{2\mu/\sigma^2} - 1}{(P^*_{\text{TP}} / P^*_{\text{SL}})^{2\mu/\sigma^2} - 1}$$

**Expected trade duration**:

$$\mathbb{E}[\tau^*] = \frac{1}{\mu}\left[\frac{p_{\text{TP}} \ln(P^*_{\text{TP}} / P_0) + p_{\text{SL}} \ln(P_0 / P^*_{\text{SL}})}{1}\right]$$

---

### 3.3 Risk-Reward Ratio

The **Risk-Reward Ratio (RR)** for a scalp trade:

$$RR = \frac{|P_{\text{TP}} - P_{\text{entry}}|}{|P_{\text{SL}} - P_{\text{entry}}|} = \frac{W}{L}$$

**Minimum required win rate** to break even given $RR$:

$$p^* = \frac{1}{1 + RR}$$

**Expected profit per trade** as a function of win rate $p$ and RR:

$$\mathbb{E}[\Pi] = p \cdot W - (1-p) \cdot L = p \cdot RR \cdot L - (1-p) \cdot L = L \bigl[p(1+RR) - 1\bigr]$$

**Profit factor** (total gross profit / total gross loss):

$$PF = \frac{p \cdot W}{(1-p) \cdot L} = \frac{p \cdot RR}{1 - p}$$

---

### 3.4 Position Sizing: Kelly Criterion

The **full Kelly fraction** of capital to risk per trade:

$$f^* = \frac{p \cdot w - q \cdot l}{w \cdot l} = \frac{p}{l} - \frac{q}{w}$$

where:
- $p = $ win probability, $q = 1 - p = $ loss probability
- $w = $ net profit per unit staked on a win
- $l = $ net loss per unit staked on a loss (positive number)

**Simplified Kelly** when $w = l$ (symmetric payoff):

$$f^* = p - q = 2p - 1$$

**Fractional Kelly** (used in practice to reduce variance):

$$f_{\text{frac}} = \kappa \cdot f^*, \qquad 0 < \kappa \leq 1$$

**Expected logarithmic growth rate** (Kelly objective):

$$G(f) = p \ln(1 + fw) + q \ln(1 - fl)$$

**Maximisation condition** $G'(f^*) = 0$ yields:

$$\frac{pw}{1 + f^*w} = \frac{ql}{1 - f^*l}$$

**Optimal position size** in units:

$$N^* = \frac{f^* \cdot \text{Account}}{\text{Price per unit} \cdot \text{Contract size}}$$

---

## 4. Execution Mathematics

### 4.1 Slippage Models

**Realised slippage** on a single trade:

$$\text{Slippage} = P_{\text{fill}} - P_{\text{decision}}$$

**Expected total execution cost decomposition**:

$$\mathbb{E}[\text{Cost}] = \underbrace{\frac{1}{2}S}_{\text{half-spread}} + \underbrace{I(Q)}_{\text{market impact}} + \underbrace{C_{\text{fees}}}_{\text{commissions/fees}}$$

**VWAP slippage** relative to volume-weighted average price:

$$\text{VWAP Slippage} = \frac{P_{\text{fill}} - \text{VWAP}}{\text{VWAP}}$$

**Implementation shortfall** (Perold 1988):

$$IS = \underbrace{P_{\text{fill}} - P_0}_{\text{price impact}} + \underbrace{P_0 - P_{\text{decision}}}_{\text{delay cost}} + \underbrace{C_{\text{fees}}}_{\text{fees}}$$

where $P_0$ is the mid-quote at the time of order submission and $P_{\text{decision}}$ is the mid-quote at decision time.

---

### 4.2 Execution Probability

**Fill probability** for a limit order at price $p$ placed at time $t$, given mid-price $M(t)$ and volatility $\sigma$:

$$P(\text{fill}\mid p, t) = P\bigl(\min_{s \in [t, t+T]} P(s) \leq p \bigr) \quad \text{(for a buy limit order)}$$

For GBM, the **first-passage probability** to level $p < M(t)$ within time $T$:

$$P(\text{fill}) = \Phi\!\left(\frac{\ln(p/M) - \mu T}{\sigma\sqrt{T}}\right) + e^{2\mu \ln(p/M)/\sigma^2} \Phi\!\left(\frac{\ln(p/M) + \mu T}{\sigma\sqrt{T}}\right)$$

where $\Phi$ is the standard normal CDF.

**Fill rate as a function of queue position** $q$ at a given price level (Gould et al. 2013):

$$P(\text{fill} \mid q) = 1 - \exp\!\left(-\frac{\lambda_{\text{trade}} \cdot T}{D(p) + q}\right)$$

where $\lambda_{\text{trade}}$ is the trade arrival rate and $D(p)$ is the depth at level $p$.

---

### 4.3 Fill Rate Calculations

**Expected fill volume** for a limit order of size $Q$ placed at depth $\delta$:

$$\mathbb{E}[\text{Fill}] = Q \cdot P(\text{fill}) = Q \cdot \left(1 - e^{-\mu_{\delta} T}\right)$$

where $\mu_\delta$ is the execution intensity at depth $\delta$.

**Partial fill probability** for order size $Q$ when available volume at level is $V$:

$$P(\text{partial fill} = v) = \frac{\binom{V}{v} \binom{N-V}{Q-v}}{\binom{N}{Q}}, \quad v \in \{0, 1, \ldots, \min(Q,V)\}$$

**Fill rate ratio** (realised fills / submitted orders):

$$\text{FR} = \frac{\text{Number of filled orders}}{\text{Number of submitted orders}}$$

---

### 4.4 Timing Optimization using Optimal Stopping Theory

**Optimal stopping problem** — find the stopping time $\tau^*$ that maximises:

$$V(x) = \sup_{\tau \geq 0} \mathbb{E}^x\!\left[e^{-r\tau} g(X_\tau)\right]$$

where $r > 0$ is the discount rate and $g(\cdot)$ is the payoff function.

**Free-boundary condition** — the optimal stopping region $\mathcal{S}^*$ satisfies:

$$\mathcal{S}^* = \{x : V(x) = g(x)\}$$

**Continuation region**: $\mathcal{C} = \{x : V(x) > g(x)\}$, where $V$ satisfies:

$$\mathcal{L} V(x) - r V(x) = 0, \quad x \in \mathcal{C}$$

**Smooth-pasting conditions** at boundary $x^*$:

$$V(x^*) = g(x^*), \qquad V'(x^*) = g'(x^*)$$

**Wald's identity** for expected gain at stopping:

$$\mathbb{E}[S_{\tau^*}] = \mathbb{E}[\tau^*] \cdot \mu$$

---

## 5. Profitability Analysis

### 5.1 Scalp Profit Function

**Gross profit** of a sequence of $N$ scalp trades:

$$\Pi_{\text{gross}} = \sum_{i=1}^{N} d_i \left(P_{\text{exit},i} - P_{\text{entry},i}\right) \cdot Q_i$$

where $d_i \in \{+1, -1\}$ denotes long/short direction and $Q_i$ is the position size.

**Net profit** after costs:

$$\Pi_{\text{net}} = \Pi_{\text{gross}} - \sum_{i=1}^{N} \left[C_{\text{spread},i} + C_{\text{commission},i} + C_{\text{slippage},i}\right]$$

**Per-trade net profit**:

$$\pi_i = d_i \left(P_{\text{exit},i} - P_{\text{entry},i}\right) Q_i - S_i Q_i - 2 f_i Q_i$$

where $S_i$ is the spread paid and $f_i$ is the commission per unit.

---

### 5.2 Expected Value of PnL

**Expected PnL per trade**:

$$\mathbb{E}[\text{PnL}] = p \cdot W - q \cdot L$$

where $p$ = win probability, $q = 1-p$ = loss probability, $W$ = average win size, $L$ = average loss size.

**Expected PnL after $N$ trades**:

$$\mathbb{E}[\Pi_N] = N \cdot \mathbb{E}[\text{PnL}] = N(pW - qL)$$

**Variance of PnL per trade**:

$$\text{Var}[\text{PnL}] = p(1-p)(W+L)^2 + p \sigma_W^2 + q \sigma_L^2$$

where $\sigma_W^2$ and $\sigma_L^2$ are the variances of win and loss sizes respectively.

**Variance of total PnL** after $N$ independent trades:

$$\text{Var}[\Pi_N] = N \cdot \text{Var}[\text{PnL}]$$

---

### 5.3 Win Rate and Average Win/Loss

**Expected value per trade** expressed using win/loss rates:

$$\mathbb{E}[\text{Trade}] = (\text{Win\%} \times \text{AvgW}) - (\text{Loss\%} \times \text{AvgL})$$

**Breakeven win rate** given average win $W$ and average loss $L$:

$$p^* = \frac{L}{W + L}$$

**Profit factor**:

$$PF = \frac{p \cdot \text{AvgW}}{(1-p) \cdot \text{AvgL}} = \frac{\text{Gross Profit}}{\text{Gross Loss}}$$

A strategy is profitable if and only if $PF > 1$.

**Edge ratio** (comparing win/loss sizes relative to adverse excursion):

$$\text{ER} = \frac{\text{AvgW} / \text{MAE}_{\text{wins}}}{\text{AvgL} / \text{MFE}_{\text{losses}}}$$

---

### 5.4 PnL Distribution and Statistical Moments

**Central moments** of the trade PnL distribution:

$$\mu_1 = \mathbb{E}[\pi] = pW - qL \quad \text{(mean)}$$

$$\mu_2 = \mathbb{E}[(\pi - \mu_1)^2] \quad \text{(variance)}$$

$$\mu_3 = \mathbb{E}[(\pi - \mu_1)^3] \quad \text{(skewness numerator)}$$

$$\mu_4 = \mathbb{E}[(\pi - \mu_1)^4] \quad \text{(kurtosis numerator)}$$

**Standardised skewness and excess kurtosis**:

$$\gamma_1 = \frac{\mu_3}{\mu_2^{3/2}}, \qquad \gamma_2 = \frac{\mu_4}{\mu_2^2} - 3$$

**Normal approximation** (by CLT) for total PnL over $N$ trades:

$$\Pi_N \approx \mathcal{N}\!\left(N \mu_1,\; N \mu_2\right)$$

**Sharpe Ratio** (annualised, assuming $n$ trades per day, $D$ trading days):

$$SR = \frac{\sqrt{nD} \cdot \mu_1}{\sqrt{\mu_2}}$$

**Sortino Ratio** (penalises only downside deviation $\sigma^-$):

$$\text{Sortino} = \frac{\mu_1 - r_f}{\sigma^-}, \qquad \sigma^- = \sqrt{\mathbb{E}\!\left[\min(\pi - r_f, 0)^2\right]}$$

---

## 6. Risk Management for Scalping

### 6.1 Drawdown Mathematics

**Drawdown** at time $t$ relative to running peak:

$$DD(t) = \frac{V(t) - \max_{0 \leq s \leq t} V(s)}{\max_{0 \leq s \leq t} V(s)}$$

**Maximum drawdown (MDD)** over horizon $[0,T]$:

$$MDD = \max_{0 \leq t \leq T} DD(t) = \frac{\max_{0 \leq s \leq T} V(s) - \min_{s \leq t \leq T} V(t)}{\max_{0 \leq s \leq T} V(s)}$$

**Expected maximum drawdown** for a GBM process with drift $\mu$ and vol $\sigma$ over $T$ periods:

$$\mathbb{E}[MDD] \approx \sigma \sqrt{T} \cdot \left(\sqrt{2 \ln n} - \frac{\ln \ln n + \ln 4\pi}{2\sqrt{2 \ln n}}\right) - \mu T$$

**Calmar Ratio** (return relative to MDD):

$$\text{Calmar} = \frac{R_{\text{annual}}}{|MDD|}$$

**Recovery factor**:

$$RF = \frac{\Pi_{\text{net}}}{|MDD|}$$

---

### 6.2 Maximum Adverse Excursion (MAE)

**MAE for a long trade** entered at $P_{\text{entry}}$:

$$MAE = \min_{t \in [\tau_{\text{entry}}, \tau_{\text{exit}}]} P(t) - P_{\text{entry}} \leq 0$$

**MAE for a short trade**:

$$MAE = P_{\text{entry}} - \max_{t \in [\tau_{\text{entry}}, \tau_{\text{exit}}]} P(t) \leq 0$$

**Maximum Favourable Excursion (MFE)** for a long trade:

$$MFE = \max_{t \in [\tau_{\text{entry}}, \tau_{\text{exit}}]} P(t) - P_{\text{entry}} \geq 0$$

**MAE-based stop optimisation** — the optimal stop $s^*$ maximises expected value:

$$s^* = \arg\max_{s \leq 0} \left\{ p(s) \cdot \mathbb{E}[W \mid \text{win}] - (1-p(s)) \cdot |s| \right\}$$

where $p(s) = P(MAE > s)$ is the probability that the trade does not hit stop $s$.

---

### 6.3 Risk per Trade

**Fixed fractional risk** per trade as a fraction $f$ of account equity $E$:

$$\text{RiskPerTrade} = f \cdot E$$

**Position size** derived from the fixed fractional method:

$$Q = \frac{\text{RiskPerTrade}}{|P_{\text{entry}} - P_{\text{SL}}| + \text{Spread}}$$

**Dollar risk** for a futures contract with tick size $\Delta$ and tick value $V_{\text{tick}}$:

$$\text{DollarRisk} = \frac{|P_{\text{entry}} - P_{\text{SL}}|}{\Delta} \cdot V_{\text{tick}} \cdot Q$$

**Maximum position size** constrained by margin requirement $m$:

$$Q_{\max} = \frac{E \cdot f_{\text{margin}}}{m \cdot P}$$

---

### 6.4 Portfolio Heat

**Total portfolio risk** (sum of individual trade risks):

$$\text{Portfolio Heat} = \sum_{i=1}^{k} \text{Risk}_i = \sum_{i=1}^{k} f_i \cdot E$$

**Correlated risk adjustment** — if trades have pairwise correlation $\rho_{ij}$:

$$\text{Portfolio Variance} = \sum_{i=1}^{k} \sigma_i^2 + 2\sum_{i < j} \rho_{ij}\, \sigma_i\, \sigma_j$$

In vector notation with covariance matrix $\Sigma$ and position vector $\mathbf{w}$:

$$\text{Portfolio Risk} = \sqrt{\mathbf{w}^\top \Sigma\, \mathbf{w}}$$

**Heat limit** — scalper caps total open risk:

$$\text{Portfolio Heat} \leq H_{\max}, \qquad H_{\max} \in [0.02E,\; 0.10E]$$

**Daily loss limit** (circuit breaker):

$$\Pi_{\text{daily}} \geq -L_{\max}, \qquad L_{\max} = \alpha \cdot E$$

If $\Pi_{\text{daily}} < -L_{\max}$, trading ceases for the day.

---

## 7. Volatility Models for Scalping

### 7.1 Intraday Volatility Estimation

**Parkinson (1980) high-low range estimator** of intraday volatility:

$$\hat{\sigma}^2_{\text{Park}} = \frac{1}{4 \ln 2} \left(\ln \frac{H}{L}\right)^2$$

where $H$ and $L$ are the intraday high and low prices.

**Garman–Klass (1980) estimator** incorporating open $O$, close $C$, high $H$, low $L$:

$$\hat{\sigma}^2_{\text{GK}} = 0.5 \left(\ln\frac{H}{L}\right)^2 - (2\ln 2 - 1)\left(\ln\frac{C}{O}\right)^2$$

**Rogers–Satchell (1991) estimator** (robust to drift):

$$\hat{\sigma}^2_{\text{RS}} = \ln\frac{H}{C}\ln\frac{H}{O} + \ln\frac{L}{C}\ln\frac{L}{O}$$

**Average intraday volatility** over $N$ periods:

$$\bar{\sigma}_{\text{intraday}} = \sqrt{\frac{1}{N}\sum_{k=1}^{N} \hat{\sigma}^2_k}$$

**Intraday volatility pattern** — the U-shaped (diurnal) effect, modelled as:

$$\sigma(t) = \bar{\sigma} \cdot f(t), \qquad \frac{1}{T}\int_0^T f(t)\,dt = 1$$

where $f(t)$ is a deterministic diurnal adjustment factor.

---

### 7.2 Volatility Clustering

**Autocorrelation of squared returns** — evidence for volatility clustering:

$$\rho_{r^2}(h) = \text{Corr}(r_t^2,\, r_{t-h}^2) > 0 \quad \text{for small lags } h$$

**Volatility persistence** — Hurst exponent $H > 0.5$ indicates long memory:

$$\mathbb{E}[|r_{t+h} - r_t|^2] \sim h^{2H}$$

**ARCH($q$) model** (Engle 1982):

$$r_t = \sigma_t \varepsilon_t, \qquad \varepsilon_t \overset{\text{i.i.d.}}{\sim} \mathcal{N}(0,1)$$

$$\sigma_t^2 = \omega + \sum_{j=1}^{q} \alpha_j r_{t-j}^2$$

---

### 7.3 GARCH(1,1) Model

**GARCH(1,1) conditional variance equation**:

$$\sigma_t^2 = \omega + \alpha\, r_{t-1}^2 + \beta\, \sigma_{t-1}^2$$

with constraints:
- $\omega > 0$ (baseline variance)
- $\alpha \geq 0$ (ARCH effect — sensitivity to past shocks)
- $\beta \geq 0$ (GARCH effect — persistence of variance)
- $\alpha + \beta < 1$ (covariance stationarity)

**Mean equation** (AR(1) for returns):

$$r_t = \mu + \phi r_{t-1} + \varepsilon_t, \qquad \varepsilon_t = \sigma_t z_t, \quad z_t \overset{\text{i.i.d.}}{\sim} \mathcal{N}(0,1)$$

**Long-run (unconditional) variance**:

$$\bar{\sigma}^2 = \frac{\omega}{1 - \alpha - \beta}$$

**Multi-step ahead variance forecast**:

$$\mathbb{E}_t[\sigma_{t+h}^2] = \bar{\sigma}^2 + (\alpha + \beta)^{h-1}\left(\sigma_{t+1}^2 - \bar{\sigma}^2\right)$$

**Log-likelihood function** for GARCH(1,1) estimation:

$$\ell(\theta) = -\frac{T}{2}\ln(2\pi) - \frac{1}{2}\sum_{t=1}^{T} \left[\ln\sigma_t^2 + \frac{r_t^2}{\sigma_t^2}\right]$$

---

### 7.4 Stochastic Volatility

**Heston (1993) stochastic volatility model**:

$$dP(t) = \mu P(t)\,dt + \sqrt{v(t)}\,P(t)\,dW_1(t)$$

$$dv(t) = \kappa\bigl(\theta - v(t)\bigr)dt + \xi\sqrt{v(t)}\,dW_2(t)$$

$$\text{Cov}(dW_1, dW_2) = \rho\,dt$$

where:
- $v(t) = \sigma^2(t)$: instantaneous variance
- $\kappa > 0$: mean-reversion speed
- $\theta > 0$: long-run variance
- $\xi > 0$: volatility of volatility
- $\rho \in (-1,1)$: leverage correlation

**Feller condition** (ensures $v(t) > 0$ a.s.):

$$2\kappa\theta > \xi^2$$

**SABR model** (Hagan et al. 2002) — commonly used for gold options:

$$dF(t) = \hat{\sigma}(t)\, F(t)^\beta\, dW_1(t)$$

$$d\hat{\sigma}(t) = \nu\, \hat{\sigma}(t)\, dW_2(t)$$

$$\text{Cov}(dW_1, dW_2) = \rho\,dt$$

**Approximate implied volatility** under SABR (Hagan formula):

$$\sigma_B(K, F) \approx \frac{\alpha}{(FK)^{(1-\beta)/2}} \cdot \frac{z}{\chi(z)} \left[1 + \left(\frac{(1-\beta)^2}{24}\frac{\alpha^2}{(FK)^{1-\beta}} + \frac{\rho\beta\alpha\nu}{4(FK)^{(1-\beta)/2}} + \frac{2-3\rho^2}{24}\nu^2\right)T\right]$$

where $z = \frac{\nu}{\alpha}(FK)^{(1-\beta)/2}\ln\!\frac{F}{K}$, $\chi(z) = \ln\!\frac{\sqrt{1 - 2\rho z + z^2} + z - \rho}{1-\rho}$, $\alpha = \hat{\sigma}(0)$ is the initial volatility level, and $\nu$ is the volatility of volatility.

---

## Summary Table of Key Scalping Formulas

| Quantity | Formula |
|---|---|
| Bid-Ask Spread | $S(t) = A(t) - B(t)$ |
| Price Dynamics | $dP = \mu P\,dt + \sigma P\,dW$ |
| Jump Intensity | $\mathbb{E}[N(T)] = \lambda T$ |
| Price Impact | $I(Q) = \alpha|Q|^\beta \cdot \text{sgn}(Q)$ |
| Kelly Fraction | $f^* = (pw - ql)/(wl)$ |
| Risk-Reward Ratio | $RR = |TP - \text{Entry}| / |SL - \text{Entry}|$ |
| Expected PnL | $\mathbb{E}[\text{PnL}] = p \cdot W - q \cdot L$ |
| Drawdown | $DD(t) = (V_{\text{peak}} - V(t)) / V_{\text{peak}}$ |
| MAE (long) | $MAE = \min P(t) - P_{\text{entry}}$ |
| Risk per Trade | $\text{RiskPerTrade} = f \cdot E$ |
| GARCH Variance | $\sigma_t^2 = \omega + \alpha r_{t-1}^2 + \beta\sigma_{t-1}^2$ |
| Intraday Volatility | $\hat{\sigma}^2 = \frac{(\ln H/L)^2}{4\ln 2}$ |
| Sharpe Ratio | $SR = \sqrt{nD}\,\mu_1 / \sqrt{\mu_2}$ |
| Implementation Shortfall | $IS = (P_{\text{fill}} - P_0) + (P_0 - P_{\text{dec}}) + C_{\text{fees}}$ |

---

*Document prepared for gold scalping strategy research. All mathematical notation follows standard financial mathematics conventions. Model parameters should be calibrated to current XAU/USD microstructure data.*
