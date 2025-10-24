---
editor_options: 
  markdown: 
    wrap: 72
---

# Revised ES sensitivity model based on [Kroiss et al](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0208336#sec027)

# Latent disease prevalence

Note, I think we could just as easily make this a model of incidence,
and then add delay distributions on case detections vs ES detections.

Model assumes that underlying prevalence of enterovirus is highly
cyclical by month. It also assumes that the pattern in cyclic behavior
will vary by location (province and district-levels).

We can write a model for virus prevalence over time in each province $i$
and district $j$, $p_{ij}(t)$ as:

$$
p_{i,j}(t) = \beta_0 + \beta_i + c_{i,j} + f_{global,t}(monthofyear) + f_{i,j,t}(monthofyear) + \alpha p_{i,j}(t-1) 
$$

Where $\beta_0$ is a global intercept (baseline prevalence of virus in
the country), $\beta_i$ is a district-level random intercept, $c_{i,j}$
is a province-level random intercept, $f_{global,t}(monthofyear)$ is the
global effect of month of the year $f_{i,j, t}(monthofyear)$ is the
province level effect of month of year and $\alpha p_{i,j}(t-1)$ is an
autoregressive term describing the dependence on the estimate at the
previous time point

Note we can modify this however we want, the point will just be we want
something that can vary, across month of year, across years, and across
provinces.

Then we have our two types of observations, both of which are functions
of the prevalence in each province at each time point.

First, for the clinical cases, we will assume that the number of cases
detected in each province and district at time $t$, $y_{i,j,t}$ among
those tested (NP-AFP individuals), $n_{i,j,t}$ follows a binomial
obervation model:

$$ 
y_{i,j,t} \sim Bin(n_{i,j,t}, p_{i,j}(t))
$$

This assumes that the prevalence of virus among NP-AFP cases is
representative of virus prevalence in the population at large, and there
aren't systematic biases in the individuals presenting with NP-AFP
cases. This differs from if we had clinical cases presenting with
symptoms, in which case, we would want to model from incidence --\>
healthcare seeking behavior.

For ES detections, we can write a model that describes the latent
probability of detection in each site $k$ in district $i$ and province
$j$ over time, $q_{i,j,k}(t)$. This will be a function of the latent
prevalence $p_{i,j}(t)$.

Following what they have written (but I do not understand so we should
revisit this)

$$
log(-log(1-q_{i,j,k}(t))) = \kappa_{0} + \kappa_{k} + f_{k,t}(monthofyear) + \epsilon_{k}log(p_{i,j}(t)) 
$$

Where $\kappa_{0}$ is some transformed baseline probability of detecting
a positive (e.g. false positive rate) $\kappa_{k}$ is some transformed
site specific intercept $f_{k,t}(monthofyear)$ is the site-specific
effect of the month of the year (this would acount for temperature,
rainfall, etc and its impact on the site's sensitivity) and \$
\epsilon*{k}log(p\*{i,j}(t))\$ is the impact of prevelance on detection
probability

Lastly, the observation model for ES detections, $D_{k,t}$, which will
be in the form of either presence or absence in a particular site $k$
as:

$$
D_{k,t} \sim Bernoulli(q_{i,j,k}(t))
$$
