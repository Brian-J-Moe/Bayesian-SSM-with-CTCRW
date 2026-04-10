# Bayesian State-Space Model with Continuous Time Correlated Random Walk
Brian J. Moe

<sup>1</sup> Fish and Wildlife Research Institute, Florida Fish and
Wildlife Conservation Commission

<sup>✉</sup> Correspondence: [Brian J. Moe
\<brian.moe@myfwc.com\>](mailto:brian.moe@myfwc.com)



## Movement model

Using a hierarchical Bayesian state-space model (BSSM) implementing a
continuous-time correlated random walk (CTCRW), we can model the latent
movement process as a function of velocity (speed and direction) through
continuous time with autocorrelation of consecutive velocity states
([Michelot et al. 2019](#ref-michelot_etal_2019); [Michelot and
Blackwell 2019](#ref-michelot_blackwell_2019)).

### Process model

The continuous-time process model for individual $i$ follows an
Orstein-Ulenbeck process with an additive drift force:

<span id="eq-process">$$
dV_i(t) = \left[-\beta_iV_i(t) + f(\textbf{x}, t)\right]dt + \sigma_idW(t)
 \qquad(1)$$</span>

where $V_i(t)$ is the velocity of individual $i$ at time $t$, $\beta_i$
is the velocity autocorrelation decay rate, $\sigma_i$ is the velocity
diffusion coefficient, $W(t)$ is the standard Weiner process, and
$f(\textbf{x}, t)$ is a position- and time-dependent drift force.
Position changes as:

<span id="eq-position">$$
dX_i(t)=V_i(t)dt
 \qquad(2)$$</span>

### Drift specification

The drift force comprises three additive components representing
behavioral thermoregulation, depth preference, and creek refuge
attraction:

<span id="eq-drift">$$
f(\textbf{x}, t) = f_{temp} + f_{depth} + f_{creek}
 \qquad(3)$$</span>

The thermal drift force was was modeled as a bidirectional restoring
force:

<span id="eq-tempf">$$
f_{temp} = \alpha\left(T_{opt} - T(\textbf{x},t) \right) \cdot \nabla T(\textbf{x}, t)
 \qquad(4)$$</span>

where $\alpha$ is the strength of the thermoregulatory response,
$T_{opt}$ is the temperature at which no thermal drift occurs, and
$\nabla T(\textbf{x}, t)$ is the spatiotemporal temperature gradient.
Temperature gradients were estimated by building temperature maps from
the means of hourly recorded temperatures at each receiver. Temperatures
at points between receivers were interpolated using the Kriging method
with the `gstat` R package ([Pebesma 2004](#ref-gstat)).

The depth drift was modeled as a restoring force with ontogenetic shift:

<span id="eq-depthf">$$
f_{depth} = \omega (D_{pref}(age_t)-D(\textbf{x}, t)) \cdot \nabla D(\textbf{x})
 \qquad(5)$$</span>

<span id="eq-depthpref">$$
D_{pref}(age_t) = d_\infty - (d_\infty - d_0) e^{-\kappa_d age_i}
 \qquad(6)$$</span>

where $\omega$ is the strength of depth-seeking behavior,
$D(\textbf{x}, t)$ is the tidally-adjusted water depth at position
$\textbf{x}$ and time $t$, $\nabla D(\textbf{x})$ is the spatial
gradient of the bathymetric surface, $D_{pref}(age_t)$ is the
age-dependent depth preference, $d_0$ is the depth preference of
neonatal individuals, $d_\infty$ is the asymptotic depth preference of
large juveniles and adults within the Glover Bight system, $\kappa_d$ is
the rate of change of ontogenetic depth preference, and $age_t$ is the
estimated age at time $t$ .

The bathymetric surface map downloaded from the National Ocean and
Atmespheric Administration Digital Coast database. We use the Nation
Center for Environemntal Information continuously updated digital
elevation model (CUDEM) which mapped the bathymetric surface at a ninth
arch-second resolution (3m) ([Cooperative Institute for Research in
Environmental Sciences (CIRES) at the University of Colorado, Boulder
2014](#ref-cudem_data); [Amante et al. 2023](#ref-cudem)). To account
for uncertainty of animal positions and the prediction error of the GAM
used to estimate animal HPEm, the bathymetric surface was smoothed to a
5m resultion using a Gaussian kernel. Spatial depth gradients were
calculated using the finite differences method with the R package
`terra` ([Hijmans et al. 2026](#ref-terra)).

### Discrete-time state transition

### Hierarchical structure

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-cudem" class="csl-entry">

Amante, C.J., Love, M.R., Carignan, K.S., Sutherland, M.G., Lim, E.,
Warnken, J., and Weisberg, R.H. 2023. Continuously updated digital
elevation models (CUDEMs) to support coastal inundation modeling. Remote
Sensing **15**(6): 1702.
doi:[10.3390/rs15061702](https://doi.org/10.3390/rs15061702).

</div>

<div id="ref-cudem_data" class="csl-entry">

Cooperative Institute for Research in Environmental Sciences (CIRES) at
the University of Colorado, Boulder. 2014. Continuously updated digital
elevation model (CUDEM) – 1/9 arc-second resolution
bathymetric-topographic tiles. NOAA National Centers for Environmental
Information.
doi:[10.25921/ds9v-ky35](https://doi.org/10.25921/ds9v-ky35).

</div>

<div id="ref-terra" class="csl-entry">

Hijmans, R.J., Brown, A., and Barbosa, M. 2026. Terra: Spatial data
analysis. Available from <https://github.com/rspatial/terra>.

</div>

<div id="ref-michelot_blackwell_2019" class="csl-entry">

Michelot, T., and Blackwell, P.G. 2019. State-switching continuous-time
correlated random walks. Methods in Ecology and Evolution **10**(5):
637–649. Wiley Online Library.

</div>

<div id="ref-michelot_etal_2019" class="csl-entry">

Michelot, T., Gloaguen, P., Blackwell, P.G., and Étienne, M.-P. 2019.
The langevin diffusion as a continuous-time model of animal movement and
habitat selection. Methods in ecology and evolution **10**(11):
1894–1907. Wiley Online Library.

</div>

<div id="ref-gstat" class="csl-entry">

Pebesma, E.J. 2004. Multivariable geostatistics in S: The gstat package.
Computers & Geosciences **30**: 683–691. Available from
<https://doi.org/10.1016/j.cageo.2004.03.012>.

</div>

</div>
