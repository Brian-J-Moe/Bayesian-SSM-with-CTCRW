# Following the Drift: A Bayesian state-space CTCRW with environmental drivers and refuge attraction




## Movement model

Using a hierarchical Bayesian state-space model (BSSM) implementing a continuous-time correlated random walk (CTCRW), we can model the latent movement process as a function of velocity (speed and direction) through continuous time ([Michelot et al. 2019](#ref-michelot_etal_2019); [Michelot and Blackwell 2019](#ref-michelot_blackwell_2019)). The CTCRW model assumes velocities over a given time interval are correlated to the that of the previous time interval. In other words, speed and directional movement tend to persist through time. The velocity process is governed by a velocity autocorrelation decay rate which defines how quickly autocorrelation decays over time and a velocity diffusion coefficient which describes the spread of velocity around the global velocity mean (which is assumed to be zero).

### Process model

The velocity model for individual $i$ follows an Orstein-Ulenbeck process with an additive drift force:

<span id="eq-process">$$
dV_i(t) = \left[-\beta_iV_i(t) + f(\textbf{x}, t)\right]dt + \sigma_idW(t)
 \qquad(1)$$</span>

where $V_i(t)$ is the velocity of individual $i$ at time $t$, $\beta_i$ is the velocity autocorrelation decay rate, $\sigma_i$ is the velocity diffusion coefficient, $W(t)$ is the standard Weiner process, and $f(\textbf{x}, t)$ is a position- and time-dependent drift force. From this, position changes as:

<span id="eq-position">$$
dX_i(t)=V_i(t)dt
 \qquad(2)$$</span>

### Drift specification

The drift force comprises three additive components representing behavioral thermoregulation, depth preference, and creek refuge attraction:

<span id="eq-drift">$$
f(\textbf{x}, t) = f_{temp} + f_{depth} + f_{creek}
 \qquad(3)$$</span>

#### Thermal restoring force

The thermal drift force was was modeled as a bidirectional restoring force:

<span id="eq-tempf">$$
f_{temp} = \alpha\left(T_{opt} - T(\textbf{x},t) \right) \cdot \nabla T(\textbf{x}, t)
 \qquad(4)$$</span>

where $\alpha$ is the strength of the thermoregulatory response, $T_{opt}$ is the temperature at which no thermal drift occurs, and $\nabla T(\textbf{x}, t)$ is the spatiotemporal temperature gradient. Temperature gradients were estimated by building temperature maps from the means of hourly recorded temperatures at each receiver. Temperatures at points between receivers were interpolated using the Kriging method with the `gstat` R package ([Gräler et al. 2016](#ref-gstat)).

#### Depth restoring force with ontogenetic shift

The depth drift was modeled as a restoring force with ontogenetic shift:

<span id="eq-depthf">$$
f_{depth} = \omega (D_{pref}(age_t)-D(\textbf{x}, t)) \cdot \nabla D(\textbf{x})
 \qquad(5)$$</span>

<span id="eq-depthpref">$$
D_{pref}(age_t) = d_\infty - (d_\infty - d_0) e^{-\kappa_d age_i}
 \qquad(6)$$</span>

where $\omega$ is the strength of depth-seeking behavior, $D(\textbf{x}, t)$ is the tidally-adjusted water depth at position $\textbf{x}$ and time $t$, $\nabla D(\textbf{x})$ is the spatial gradient of the bathymetric surface, $D_{pref}(age_t)$ is the age-dependent depth preference, $d_0$ is the depth preference of neonatal individuals, $d_\infty$ is the asymptotic depth preference of large juveniles and adults within the Glover Bight system, $\kappa_d$ is the rate of change of ontogenetic depth preference, and $age_t$ is the estimated age at time $t$ .

The bathymetric surface map downloaded from the National Ocean and Atmospheric Administration Digital Coast database. We use the Nation Center for Environemntal Information continuously updated digital elevation model (CUDEM) which mapped the bathymetric surface at a ninth arch-second resolution (3m) ([CIRES 2014](#ref-cudem_data); [Amante et al. 2023](#ref-cudem)). To account for uncertainty of animal positions and the prediction error of the GAM used to estimate animal HPEm, the bathymetric surface was smoothed to a 5m resolution using a Gaussian kernel. Spatial depth gradients were calculated using the finite differences method with the R package `terra` ([Hijmans et al. 2026](#ref-terra)).

#### Creek refuge attraction with size-dependent decay

Glover Bight Creek acts as a refuge from bull shark predation for yoy and small individuals. Attraction to the creek refuge was modeled as a length-dependent logistic decay:

$$
f_{creek} = \psi \cdot \omega(L_t) \cdot \nabla\tilde{C}(\text{x}),
$$ {#eq_creekf}

$$
\omega (L_t) = \frac{1}{1 + \text{exp}[\kappa_c(L_t - L_{50})]}
$$ {#eq_creekLogis}

where $\psi$ is the strength of creek attraction, $\nabla\tilde{C}(\text{x})$ is the gradient of a softplus-clamped signed distance field pointing toward the creek interior, $\omega(L_t)$ is a logistic weight which decays with length $L_t$ at time $t$, $L_{50}$ is the length at which creek affinity halves (the inflection of the curve), and $\kappa_c$ describes the steepness of the transition. The signed distance field was clamped using a softplus transformation $\tilde{C} = k^{-1}\text{ln}(1+e^{kC})$ to prevent inward drift for animals already inside the creek.

### Discrete-time state transition

### Hierarchical structure

To account for individual variation, we modeled the individual specific parameters $\beta_i$ and $\sigma_i$ hierarchically on the log scale:

<span id="eq-global_beta">$$
\text{ln}\beta_i = \text{ln}\mu_\beta + \tau_\beta z_i^{(\beta)}, \quad z_i^{(\beta)} \sim \mathcal{N}(0, 1),
 \qquad(7)$$</span>

<span id="eq-global_sigma">$$
\text{ln}\sigma_i = \text{ln}\mu_\sigma + \tau_\sigma z_i^{(\sigma)}, \quad z_i^{(\sigma)} \sim \mathcal{N}(0, 1).
 \qquad(8)$$</span>

Where $\mu_\beta$ and $\mu_\sigma$ are population level medians of $\beta_i$ and $\sigma_i$, respectively, and $\tau_\beta$ and $\tau_\sigma$ are among-individual standard deviations on the log scale.

The population level hyperparameters were given the priors

$$
\text{ln}\mu_\beta \sim \mathcal{N}(\text{ln}\widetilde{\beta}, s_\beta),
$$

$$
\text{ln}\mu_\sigma \sim \mathcal{N}(\text{ln}\widetilde{\sigma}, s_\sigma),
$$ $$
\tau_\beta \sim \text{Exponential}(\lambda_\tau), \quad \tau_\sigma \sim \text{Exponential}(\lambda_\tau).
$$

Where $\widetilde{\beta}$ and $\widetilde{\sigma}$ are the specified prior medians, $s_\beta$ and $s_\sigma$ are the prior standard deviations, and $\lambda_tau$ is the prior exponential rate. Continuous tracks were partitioned into independent segments based on intervals longer than 30 minutes. All segments of a given individual share that individuals movement parameters, allowing repeated segments to inform the same hierarchical random effects.

Models were build in Stan ([Stan developement Team 2025](#ref-stan)) using the R package `cmdstanr` ([Gabry et al. 2025](#ref-cmdstanr)).

# References

<div id="refs" class="references csl-bib-body hanging-indent" entry-spacing="0">

<div id="ref-cudem" class="csl-entry">

Amante, C.J., Love, M.R., Carignan, K.S., Sutherland, M.G., Lim, E., Warnken, J., and Weisberg, R.H. 2023. Continuously updated digital elevation models (CUDEMs) to support coastal inundation modeling. Remote Sensing **15**(6): 1702. doi:[10.3390/rs15061702](https://doi.org/10.3390/rs15061702).

</div>

<div id="ref-cudem_data" class="csl-entry">

CIRES. 2014. Continuously updated digital elevation model (CUDEM) – 1/9 arc-second resolution bathymetric-topographic tiles. NOAA National Centers for Environmental Information. doi:[10.25921/ds9v-ky35](https://doi.org/10.25921/ds9v-ky35).

</div>

<div id="ref-cmdstanr" class="csl-entry">

Gabry, J., Češnovar, R., Johnson, A., and Bronder, S. 2025. <span class="nocase">cmdstanr</span>: R interface to ’CmdStan’. Available from <https://github.com/stan-dev/cmdstanr>.

</div>

<div id="ref-gstat" class="csl-entry">

Gräler, B., Pebesma, E., and Heuvelink, G. 2016. Spatio-temporal interpolation using <span class="nocase">gstat</span>. The R Journal **8**: 204–218. Available from <https://journal.r-project.org/articles/RJ-2016-014/index.html>.

</div>

<div id="ref-terra" class="csl-entry">

Hijmans, R.J., Brown, A., and Barbosa, M. 2026. <span class="nocase">terra</span>: Spatial data analysis. doi:[10.32614/CRAN.package.terra](https://doi.org/10.32614/CRAN.package.terra).

</div>

<div id="ref-michelot_blackwell_2019" class="csl-entry">

Michelot, T., and Blackwell, P.G. 2019. State-switching continuous-time correlated random walks. Methods in Ecology and Evolution **10**(5): 637–649. doi:[10.1111/2041-210X.13275](https://doi.org/10.1111/2041-210X.13275).

</div>

<div id="ref-michelot_etal_2019" class="csl-entry">

Michelot, T., Gloaguen, P., Blackwell, P.G., and Étienne, M.-P. 2019. The Langevin diffusion as a continuous-time model of animal movement and habitat selection. Methods in ecology and evolution **10**(11): 1894–1907. doi:[10.1111/2041-210X.13154](https://doi.org/10.1111/2041-210X.13154).

</div>

<div id="ref-stan" class="csl-entry">

Stan developement Team. 2025. Stan modeling language users guide and reference manual. Available from <https://mc-stan.org/docs/2_38/>.

</div>

</div>
