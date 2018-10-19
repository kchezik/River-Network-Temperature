% Title
% ^1,2^Kyle A. Chezik; ^1^Jonathan W. Moore
% ^1^Earth to Ocean Research Group -- Simon Fraser University 8888 University Dr. Burnaby BC, Canada V5A1S6 -- ^2^778.782.9427 -- kchezik@sfu.ca

#Abstract

#Introduction

In the freshwater environment, temperature is a powerful and dynamic force that reflects complex climate and landscape interactions [@Webb:2008] and controls ecological processes [@Angilletta:2009]. As a filtered manifestation of climate, freshwater temperatures and their impacts on biota are of increasing concern with rising global temperatures. The contribution of climate change to shifting freshwater temperatures is complexified by landscape dynamics. For example, lakes are known to act as heat sinks, warming and stabilizing the thermal regimes of downstream rivers [@Wetzel:2001], while glaciers and steep slopes often contribute to stream cooling through the contribution of melt- and groundwater [@Cadbury:2008; @Lisi:2015]. The complex mixture of thermal contributions to freshwater temperature is a hierarchical aggregate of the contributing watershed (i.e., upslope area) [@Vannote:1980; @Isaak:2014]. As catchment area increases, a greater variety of landscape dynamics contribute to water's thermal profile and their relative importance shift in magnitude. This results in a heterogeneous structuring of temperature across the river network that organizes the thermal habitat and freshwater ecosystem [@Brown:2008]. By integrating over space and time, river networks exhibit stabilizing properties, dampening variation in temperature [@Moore:2015; @Steel:2016] and buffering against the changing climate [@Chezik:2017]. As a result, thermal extremes and the associated stress on organisms will vary depending on network location and the contributing landscape features.

Organisms moving through the network experience shifting thermal habitat of varied extremes and volatility [e.g., @Steel:2016]. Depending on an organisms thermal limits, network accessibility may expand and contract seasonally with lethal temperatures precluding access and sub-lethal temperatures limiting residence time. Passage through sub-lethal temperatures may be necessary to access optimal habitat and different routes through the network will confer different levels of thermal risk. Renowned for their extraordinary migrations from the Pacific Ocean to the interior of Western North America, adult Pacific Salmon (*Oncorhynchus spp.*) forge through extreme climates and landscapes in order to reproduce in their natal streams [@Quinn:2011]. A coldwater species, adult salmon are often beyond their thermal optimum during their migration leading to thermal stress that can result in pre-spawn mortality [e.g., @Hinch:2012]. Spawning throughout river systems west of the Rockie Mountain range at different times throughout the year, species and populations experience varied levels of thermal risk depending on their migration timing and route. If water temperatures increase as expected under climate change, adult migration may become increasingly strenuous possibly leading to abundance declines [@Rand:2006; @Hague:2011].

The Thompson River watershed in central British Columbia Canada is a tributary to the greater Fraser River, a system which has exhibited substantial warming over the last half century [@Morrison:2002; @Ferrari:2007]. Home to pink (*O. gorbuscha*), chinook (*O. tshawytscha*), coho (*O. kisutch*) and sockeye (*O. nerka*) salmon as well as steelhead (*O. mykiss*), the Thompson has exhibited population declines across species due to degraded ocean conditions and freshwater habitat [@labelle:2009]. Notably, the Adams River, a tributary to the southern branch of the Thompson River, hosts one of the largest runs of sockeye salmon in the world. Unfortunately, sockeye migrating to the South Thompson have exhibited an increase in pre-spawn mortality [@Young:2006; @Hinch:2012]. Linked to a shift towards earlier adult migration, energy depletion has been cited as a cause of decreased survival with warmer migration temperatures speculated to be a fundamental contributor. Similarly, increasing temperatures is believed to be contributing to the decline of coho, steelhead and chinook abundances [@Bradford:2000]. Although direct links of temperature on these species is to date unclear, sub-lethal effects have been shown to exacerbate the impacts of disease and parasites [@Gilhousen:1992] on energy density and curtail aerobic scope [@Eliason:2011], making migration more strenuous and survival less likely [@Martins:2011]. 

Understanding the various contributors to stream temperature during migration should inform routes and populations at greater risk of thermal stress and contribute a new dimension of consideration when prioritizing conservation. Using temperature data collected throughout the Thompson River watershed over four consecutive years, we built a stream network model to link climate and landscape contributors to water column temperature. Utilizing Bayesian and parametric bootstrapping methods we uniquely aggregate the probability of extreme temperatures upstream and map the thermal risk of all potential migration routes. This study reveals consistent modification of climate by landscape features that result in varied directional contributions to cumulative migratory thermal stress. Importantly, this work demonstrates the impact of limited inter-annual climatic variation on thermal stress, with implications for increased thermal risk under climate change.

#Methods

###Data

*Thompson River Watershed*

To capture the thermal heterogeneity of the Thompson River watershed, Hobo Pendant data loggers were deployed at 103 locations with the goal of capturing confluence interactions and specific landscape effects (Fig. $\ref{fig:1}$). By placing loggers in triads, upstream-stream and downstream of major confluences, we aimed to capture the relative contribution of each stream to the subsequent downstream temperatures [i.e., @Marsha:2018]. We focused logger deployment in the North Thompson because the watershed contains more distinct regions where the contributions of specific landscape features may be more pronounced. In its northern reaches, the North Thompson watershed is dominated by steep slopes, snow capped mountains and small glaciers. The West side of the watershed is perched on a high plateau, is relatively flat and dominated by lakes connected by slow moving streams. Between these two regions is a large inactive, volcanic field, largely captured within Wells Grey Provincial Park, and distinguished by canyons, rapid flows and dramatic waterfalls. The southern portion of the North Thompson is narrow, largely composed of the mainstem with many small creeks and streams flowing out of increasingly arid hills. The South Thompson watershed is dominated by lakes, which are fed by small streams tumbling out of the Canadian Rocky Mountains. These two watersheds merge in the southwest before continuing to the mainstem of the Fraser River and ultimately to the Pacific Ocean. The region is dominated by conifers, is actively logged and frequently experiences wildfires. Agriculture and ranching occur in the watershed but the topology of the region limit these industries largely to the floodplain along the mainstem with the exception of some more extensive orchard agriculture in the South Thompson watershed's Okanagan Valley.

\begin{figure}[h]
\centering
\includegraphics[width=4in]{./Water_Temperature/images/sites_pop.png}
	\caption{Stream temperature monitoring locations (blue) and chinook salmon populations (yellow) within the Thompson River watershed. The inset indicates the North (green) and South Thompson's (purple) watershed locations within the Fraser River watershed (dark grey) in British Columbia Canada (light grey). The shape of the chinook salmon population locations indicate the timing of migration through the Albion Test Fishery in early (June), mid (July) and late (August) summer, with the exception of Louis which is considered to migrate in the spring (March-May) but was grouped into the early summer period for this study (Parken et al. 2008).}
\label{fig:1}
\end{figure}

*Stream Temperature Data Collection*

Collected at two hour intervals, temperature data were gathered from July of 2014 to August of 2018, with routine annual data collection and redeployment. Annual collection was necessary to ensure the loggers remained in the water column, depleted batteries were replaced and loggers lost in the previous season were redeployed. Loggers were affixed to the landscape by attaching cable to boulders, bridge pylons or large trees expected to remain constant and unyielding to extreme river flows. The cable was then attached to a white PVC shield [e.g., @Isaak:2011], which protected loggers from debris and solar radiation that may bias temperature readings. To prevent de-watering events, loggers were weighted down by rocks and placed in deep pools protected from the direct influence of the rivers thalweg when possible. The resulting temperature data were cleaned via manual inspection, removing data believed to be associated with de-watering events that resulted in air temperature readings [e.g., @Sowder:2012]. Concerned with the most likely thermal experience, we averaged raw temperature readings to the daily level.

###Statistical Models

Interested in how temperatures change during adult salmon migration we require linked spatial and temporal models that describe thermal conditions in any given place on the river network, at any given time. To do this we first characterize temperature's temporal pattern and then subsequently model those descriptive coefficients given climatic and landscape covariates in a spatial statistical stream network model (SSNM).

*Temporal Stream Temperature Model*

To characterize the seasonal dynamics of stream temperatures we used two cyclical models in linear combination, 

\begin{linenomath*}
\begin{equation}
	y_{s,t} = Annual_{s,yr} + Season_{s,yr} + \eta_{s}, \quad
	\eta_{s} \sim \mathrm{N}(0, \sigma_{s}) \label{eq1}
\end{equation}
\end{linenomath*}

where mean daily ($t$) temperatures ($y$) at each site ($s$) are described primarily by broad seasonal dynamics ($Annual$) and modified by local seasonal conditions ($Season$). The remaining variance ($\sigma$) not captured by these models is assumed *i.i.d*.

Dominated by seasonal shifts in solar radiation due to the tilt of the earth's axis and its orbit around the sun, we describe the primary periodic dynamics of stream temperature using a cyclical model [e.g., @Shumway:2000], 

\begin{linenomath*}
\begin{equation}
	Annual_{s,yr} = \alpha_{s,yr} + A_{s,yr}\cos(2\pi\omega + \tau_{1,s}\pi)  \label{eq2}
\end{equation}
\end{linenomath*}

where the annual oscillation of temperature is described by a modified cosine curve. Shifing the cosine curve vertically captures the mean annual water temperature ($\alpha$) while amplifying the curve captures the annual range of temperatures around the mean ($A$). The frequency of the temperature cycle is captured in ($\omega$) and described by,

\begin{linenomath*}
\begin{equation}
	\omega = d_{1:\gamma}/\gamma, \label{eq3}
\end{equation}
\end{linenomath*}

where $\gamma$ is the number of observations per year and $d_{t}$ the integer location of observation $t$ in the year. Because we averaged our data to the daily level, we used the day of year as a common value of $d_{t}$ across sites, and set $\gamma$ to 365 (or 366 in leap years). Although the solstices should vaguely coincide with the extremes in annual temperature, shifts due to local climatic variation are allowed in $\tau_{1}$. Because $d_{1}$ is equivlent to January 1^st^, $\tau_{1}$ is approximately described by a value of 1 in the Northern Hemisphere or a shift in the cosine curve of $\pi$ such that the cosine curve begins near its lowest point.

Hysteresis is a common feature of annual stream temperatures, where rising spring and declining fall temperatures are not symetrical around summer's peak [@Letcher:2016], rendering a simple cosine curve a poor discriptive model on its own. Divergences from the expected pattern of temperatrure given changing levels of solar radiation are typcially a result of seasonal precipitation, where higher stream flows depress stream temperatures [@Lisi:2015]. In the northern hemisphere this is commonly observed in the spring when snow melt drives a basin wide freshet. To account for seasonal patterns in stream temperature driven by local climate, we include a double cycle cosine curve that acts to modulate the annual curve (eq. $\ref{eq2}$), 

\begin{linenomath*}
\begin{equation}
	Season_{s,yr} = \phi*\cos(2\pi\omega\tfrac{1}{2} + \tau_{2,s}\pi)  \label{eq4}
\end{equation}
\end{linenomath*}

where $\phi$ is an expansion factor that alters the strength of the hysteresis and just as in equation $\ref{eq2}$, $\tau_{2}$ shifts the curve to align with local seasonal effects. In the Thompson River watershed, $\tau_{2}$ typically corresponds with a value of 1.5 which results in warmer winter and summer stream temperatures and cooler spring and fall temperatures. This is because low flows in the summer and winter lead to a greater influence of solar radiation and ground water respectively, while snow-melt and rain dominate in the spring and fall.

Implementing this model in the probabilistic programming language Stan [@Stan:2017], allowed us to  explore the entire probability space of these parameters and propogate that variability to the subsequent spatial model. Furthermore, using a Bayesian framework allowed us to provide reasonable prior estimates and limit the parameter space over which the model must search. For instance, $\alpha$ and $A$ parameters are not reasonably expected to go below 0$\text{\textdegree}$C, therefore we limited those parameters to above zero and provided weak lognormal priors centered on 2 with a variance of 1. Similarly, given the cyclical nature of our model, $\tau$ values can take on reasonable estimates at regular intervals *ad infinitum*. Running a multi-chain (MCMC) analysis requires the chains agree on an estimate which can be difficult when multiple values lead to the same fit. To avoid disagreement we limited $\tau$ to positive values and provided weakly normal priors with a variance of 0.25 centered on 0.9 for $\tau_{1}$  and 1.5 for $\tau_{2}$. The expansion factor $\phi$ was limited to between 0 and 3, maximizing the hysteresis adjustment to no more than 3$\text{\textdegree}$C in either direction thereby limiting ths models flexibility and influence. After fitting equation $\ref{eq1}$ to the data, we visually inspected the fits at each site and year and consulted diagnostic statistics accepting Rhat value less than 1.1 and effective sampling size ratios of greater than 0.001 [@Gelman:2013]. We excluded sites and years that had poor visual fits do to a lack of data describing the critical seasonal extremes even when other diagnostics were not problematic.

*Spatial Stream Network Model*

Using cutting edge, spatial statistical stream network models [@Ver:2010], we link climate and landscape attributes to parameter estimates in equation $\ref{eq1}$ and leverage flow connected auto-covariance within and across networks to accurately predict parameters anywhere in the Thompson River watershed [@Isaak:2014]. Currently, multi-variate methods that link predictors to a variety of responses simultaneously are not available within the stream-network framework. Therefore, we built individual models for each parameter resulting in five stream network models with the general construction,

\begin{linenomath*}
\begin{equation}
\mathrm{std\_coef}_{s,yr} = \mathbf{X}\pmb{\beta} + \epsilon_{s,yr}, \quad
\epsilon \sim \mathrm{z}_{d} + \mathrm{z}_{yr} + \mathrm{z}_{s} + \mathrm{z}_{nug} \label{eq5},
\end{equation}
\end{linenomath*}

where standardized equation $\ref{eq1}$ coefficients ($\mathrm{std\_coef}$, see Table $\ref{tbl1}$) at site $s$ in year $yr$ are predicted by a matrix of climate and landscape variables $\mathbf{X}$, where the relationship between the temporal model's coefficients and these predictors are described by a vector of $\pmb{\beta}$ coefficients. The error term is decomposed into random effects ($z$) accounting for basin wide year effects ($z_{yr}$), site effects ($z_{s}$), exponentially weighted flow connected autocovariance ($z_{d}$) and residual *i.i.d* error ($z_{nug}$) [@Ver:2010].

In order to calculate predictor variables at each site, we first built a landscape network describing the branching architecture of our streams [@Theobald:2006] and delineated reach contributing areas (RCA) for each downstream confluence in the network [@Peterson:2014]. Georeferenced stream data were gathered from the British Columbia Data Catalogue made available by the Ministry of Citizens Services. Streams were subset to greater than third order, cleaned of network braiding, complex confluences, pseudo-nodes and all vectors were directed towards the outlet. Sites were snapped to the network (Fig. $\ref{fig:1}$) and center-points of each stream segment were used as prediction locations. Reach contributing areas were delineated using a 25m digital elevation model (DEM). Prior to delineation, the DEM was cleaned of topological errors and the known stream network was burned into the DEM at a 5 metre depth. All predictor variables were summarized by RCA using the *'zonal_stats'* tool in the Python package *'rasterstats'*.

Climate and landscape predictor variables known to significantly contribute to stream temperatures [e.g., @Isaak:2010], were gathered from multiple sources. Climate data were calculated using the open source tool ClimateBC [@Wang:2016], which extracts and downscales PRISM climate normal data and extrapolates to any location within British Columbia [@Daly:2008]. Sampling at 1km resolution across the Thompson River watershed, we used ClimateBC to estimate the mean annual air temperature in each year between 2014 and 2017, and calculated the amplitude around the mean as the difference between the average maximum summer temperature (June-Aug.) and the average minimum winter temperature (Dec.-Feb.). Static landscape variables included glacial coverage gathered from the Rudolph Glacier Inventory (v5.0) [@Arendt:2015], lake area calculated from polygons held in the BC Freshwater Atlas and elevation extracted from the aforementioned DEM. To characterize dominant changes to the landscape in recent decades, we extracted high resolution forest change data from the National Forest Information System. Derived from Landsat images, these data identify logging and wildfire for each year between 1985 and 2010 at a 30m resolution [@White:2017]. Climate and elevation data were averaged by RCA whereas landscape variables were summarized by total km^2^ in each RCA.

Ultimately, we are interested in the contribution of each predictor variable from the entire watershed to each observation and prediction site. Using the spatial relationships defined in our landscape network, we calculated the upslope area average of each predictor for each site. To characterize the contribution of glaciers, lakes and forest-change relative to the size of each sites catchment, we divided these variables by their catchment area, therby calculating their proportional area coverage. Visual inspection suggested a logit transformation of these percentages in order to spread the data more evenly and limit the influence of extreme data points. Climate, elevation and catchment area were all centered and scaled in order to limit the influence of the y-intercept and allow effect size comparison across variables ($\ref{tbl1}$).

\begin{table}[h!]
\centering
\begin{tabular}{ c|c|c|c }
	& Variable & Mean & SD\\
	\hline
	Response & ln($\alpha$) & & 0.29\\
	& ln($\textit{A}$) & & 0.26\\
	& $\phi$ & & 0.72\\
	& $\tau_{1}$ & & 0.04\\
	& $\tau_{2}$ & & 0.17\\
	& $\sigma$ & & 0.20\\
	\hline
	Predictor & Air Amplitude $H_{2}O$ & 29.6 & 3.55\\
	& Air Mean $H_{2}O$ & 3.11 & 1.94\\
	& ln(CA) $H_{2}O$ & 2.75 & 2.21\\
	& Elevation $H_{2}O$ (m) & 1431 & 349\\
\end{tabular}
	\caption{Mean and standard deviation (SD) values of predictor and response variables. Only the standard deviation value is reported for the response variable as these values were only standardized but not scaled. Moreover, these SD values are averages from the 250 samples drawn from equation $\ref{eq1}$ parameter posteriors. Predictor variables were both centered and scaled and the respective mean and SD values were calculated across all observed and predicted sites. Transformations of some variables by the natural log (ln) was done prior to calculating the mean and SD. $H_{2}O$ indicates predictors that constitute the average value among values averaging over the contributing area. For instance, the average contributing area elevation among predicted and observed sites was 1431 which is similar to, but not the same as, the Thompson's overal mean elevation of 1264 metres as sites are nested.}
	\label{tbl1}
\end{table}

Not all predictor variables display strong relationships with each temporal coefficient, therefore we limited predictors and interactions in each model to those that resulted in lower Akaike Information Criterion (AIC), significant gains in explained variance and improved leave one out cross validation results. Once SSNM models were constructed for each parameter, we iteratively sampled coefficient estimates from the posteriors of equation $\ref{eq1}$, fit our SSNM models to those estimates and predicted temporal coefficients at our prediction sites. This process resulted in 250 estimates of each coefficient for each year (2014-2017) at 4376 locations throught the Thompson River watershed.

To approximate the thermal risk of migrating salmon we need to define when salmon might arrive in the Thompson, how quickly they might move through the network and what temperatures are considered stressful. Migration timing of Chinook populations throughout the Thompson River watershed were estimated in 2000 and 2001 by @Parken:2008 using coded wire tags and genetic sampling. Measured at the Albion Test Fishery, appoximately 50km upstream from the Fraser River mouth, populations were grouped into three primary migration windows, early (June), mid (July) and late (August) summer. Using an estimated migration rate of 36km per day [@Salinger:2006], we calculated the number of days (~13 days) required to travel the 447km between the Albion Test Fishery and the confluence of the North and South Thopmson Rivers. Although we do not account for temperature driven changes in migration speed [e.g., @Salinger:2006], we do sample migration dates uniformly across the migration windows which should consume individual variation and generalize our findings more widely. Moreover, although we are using migration timing data for chinook salmon, we note that these summer windows are widely considered  to coincide with summer fisheries targeting sockeye and were chosen by @Parken:2008 to correspond with contemporary conservation requirements for salmon more broadly [@Bailey:2001].

Upon estimating when salmon are expected to arrive in the Thompson River system, we used our SSNM model estimates of equation $\ref{eq1}$ parameters to calculate expected frequencies of predicted daily temperatures at each prediction site in the network. Dividing the upstream distance of any point by the migration speed and adding the result to the date of arrival at the outlet, gives the expected $d$ value for each prediction site in the network. Using these values and all site predicted coeficient estimates we calculated the full distribution of mean temperature estimates at each prediction point (e.g. Fig. $\ref{fig:2}$). Considering the probability of stressful temperatures, given the mean and standard deviation estimates assuming normality, we calculated the median thermal exposure probability at each site. Using the flow connected structure of the river network, we summed all above threshold mean value (i.e., >50\% likely) occurances calculated from each site to the outlet and determined the cummulative migratory thermal exposure opportunity as the median accumulation value at each site during each migration period. Dividing the cumulative migratory thermal exposure opportunity by the total number of temperature estimates (*n* = sites $\cdot$ 250), we calculate the cummulative migratory thermal exposure probability. While the exposure opportunity strictly increases during migration, the exposure probabilty can both increase and decrease, therby depicting changing temporal and landscape changes in temperature during migration.

A key flexibility to this approach, is that it accomodates calculating thermal exposure at a variety of thresholds and spatial extents as determined by species thermal limits and behaviours. Here we use a lower theshold of 19$\text{\textdegree}$C to indicate significant but managable stress and a value of 22$\text{\textdegree}$C to indicate temperatures blocking migration per temperatures identified for sockeye, steelhead and chinook [@Richter:2005]. Ultimately, this process produces thermal exposure probabilities by migration group, in any year and at any temperature threshold of interest, at any point in the network. As such we calculated the cummulative migratory thermal exposure opportunity and spawn site thermal exposure probability for the 16 populations of chinook salmon identified by @Parken:2008 in the Thompson River watershed (Fig. $\ref{fig:1}$).

#Results

*Temporal Stream Temperature Model*

Of the 408 unique site-year combinations, we selected 290 that displayed strong fits, capturing the maximum temperatures and hysteresis present in the data (e.g., Fig. $\ref{fig:2}$). Of those the data were spread fairly evenly between years with 61 estimates in 2014, 91 in 2015, 67 in 2016 and 71 in 2017. Many sites displayed poor fits when temperature data near the annual extremes were limited and were eliminated from the spatial analysis. The selected sites in Fig. $\ref{fig:2}$ exhibit many of the characteristics that either constituted acceptance or rejection. For instance, the Barri$\`e$r River lower in the North Thompson River watershed displayed four quality years of data and acceptable fits. Meanwhile, Bridge Creek and the Shuswap River each exhibited one year, 2014 and 2016 respectively, that did not capture enough data to constitute convincing estimates of all parameters. Overal fits were strong with no apparent bias seasonally except for when peak summer temperatures were extremely 'sharp' but our estimates of variance largely captured these values at 95% confidence (e.g., Bridge Creek, 2015 Fig. $\ref{fig:2}$).

\begin{figure}[h]
\centering
\includegraphics[width=6.5in]{./Water_Temperature/images/temporal_fit.png}
	\caption{Examples of model fits to equation $\ref{eq1}$ at 4 sites in Thompson River watershed. The red ribbon describes the 95\% credibility interval of the mean temperature estimate, while the orange dashed-line describes the average 95\% variance window above the mean. The Barri$\`e$re River and Bridge Creek are tributaries to the North Thompson River (NT), while the Shuswap River drains Mabel Lake and flows into the South Thompson River (ST).}
\label{fig:2}
\end{figure}

*Spatial Stream Network Models*

Stream network (SSNM) model selection resulted in primarly strong fits that accounted for much of the variance in the temporal model paramters (eq. \ref{eq1}). Leave-one-out cross validation demonstrates the predictive capability of our models with largely strong 1:1 relationships between observed and predicted parameter values (Fig. $\ref{fig:3}$). Variance ($\sigma$) and timing ($\tau$) coefficients exhibited particularly strong predictive relationships with limited exception, likely due to their minimal variance among sites and shared values across years (i.e., eq. $\ref{eq2}$). The only covariate included in the these three SSNM models was elevation, with higher mean catchment elevations leading to lower values across parameters. In other words, higher elevation resulted in less thermal variation, a later onset of spring and delayed seasonal effects of snow and rain. The primary difference between these three models is the relative contribution of the covariate versus the correlation structures. Elevation dominated in the $\sigma$ model (46\%), while only capturing 10\% of the variation in the $\tau$ models. Site level correlation ($z_{s}$) accounted for much of the remaining variance in all three models but $\tau_{1}$ demonstrated a particularly strong effect of site (62\%) while variance in the $\tau_{2}$ parameter was primarily captured through autocorrelation ($z_{d}$) along the flow connected network (52\%). All three of these SSNM models had very little unaccounted for residual variance ($\ll$ 1\%) (Fig. $\ref{fig:3}$).

\begin{figure}[h]
\centering
\includegraphics[width=6.5in]{./Water_Temperature/images/cross_val.png}
	\caption{Evaluation of SSNM model predictive power using leave-one-out cross validation for each temporal model parameter. A 1:1 line describes a perfect fit between observed and SSNM predicted values. Parameter estimates were randomly drawn from equation \ref{eq1} parameter posteriors and predicted by SSNM models (see Fig. $\ref{fig:4}$ for covariates).}
\label{fig:3}
\end{figure}

Covariates that contributed to the efficacy of $\alpha$, $A$ and $\phi$ SSNM models varied with the exception of mean annual air temperature which was the dominant predictor across models (Fig. $\ref{fig:4}$). Water's mean annual temperature ($\alpha$) was additionally positively affected by the annual range of temperatures averaged across the contributing area (Air Amplitude), the percent catchment burned (WildFire), the catchment's mean elevation (Elevation) and the percentage of lake area (Lakes). The interaction of elevation and lake predictors was notably the only cooling influence on $\alpha$ while counter-intuitively elevation alone had a slight warming effect. Overall, these covariates captured on average 37\% of the variation in the observed data while auto-correlation along the flow connected network captured 29\%. Much of the remainder was accounted for by correlation within years (25\%) and correlation within sites (4\%) with 5\% of the variance remaining unaccounted.

\begin{figure}[h]
\centering
\includegraphics[width=6.5in]{./Water_Temperature/images/SSN.png}
	\caption{Coefficient estimates for mean annual water temperature ($\alpha$), water temperature amplitude ($A$) and seasonal hysteresis ($\phi$) SSNMs by column respectively. Top row coefficients have been scaled while the bottom row have been logit transformed making parameter estimate effect sizes comparable by row within models.}
\label{fig:4}
\end{figure}

The SSNM describing the amplitude (*A*) of water temperatures around the mean showed postive relationships with the contributing basins mean air amplitude, the catchment area and its percentage lake area. The percent catchment area covered by glaciers had a negative affect, reducing the range of annual temperatures around the mean (center, Fig. $\ref{fig:4}$). The covariates in this SSNM captured the largest amount of variance in the data (46\%) compared with covariates in all other SSNMs. The remainder was largly evenly distributed between the flow connected correlation structure ($z_{d}$ = 19\%) and correlation within sites ($z_{s}$ = 18\%).

The $\phi$ parameter controlling the strength of seasonal hysteresis was notably under-estimated by the SSNM at high values and over-estimated at low values, resulting in relatively even attribution of this parameter across sites (Fig. $\ref{fig:3}$). Neither the covariates (10\%) nor the autocorrelation structures ($z_{d}$,$z_{yr}$,$z_{s}$) dominated in accounting for variance in this parameter (20,25 and 7\% respectively). The little variance captured by covariates included negative effects of mean annual air temperatures, mean catchment elevation and percent catchment area covered by glaciers (Fig. $\ref{fig:4}$). In other words, seasonal hysteresis is reduced in colder regions, a characteristic typical of higher elevations, and further dampened by the presenece of glaciers at high elevations.

*Projected Network Thermal Exposure*

Spatial projections of thermal exposure probabilities highlight the thermal challenge of the South Thompson River's mainstem and the North Thompson River's western high plateau (Figs. $\ref{fig:5}$, S1). The potential for stressful conditions above 19$\text{\textdegree}$C were relatively elevated during the mid-summer migration and least likely during the late-summer migration (Fig. S1). Furthermore, the likelihood of stressful conditions declined between 2014 and 2017 accros migratory periods. Thermal exposure at 22$\text{\textdegree}$C was extremely unlikely across migrations and years (Fig. S2). Of those regions that had elevated stress potentials above 19$\text{\textdegree}$C, low elevation lakes constituted the most common landscape characteristic.

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{./Water_Temperature/images/maps_2015.png}
	\caption{Thermal exposure above 19 and 22\text{\textdegree}C in the North and South Thompson River watersheds during the mid-summer migration in 2015. (\textbf{left}) The median probability of experiencing above threshold temperatures throughout the network. (\textbf{center}) The median cummulative number of likely (>50\% probability) above threshold temperature estimates during upstream migration. The color bar ranges from 0 (low) to \textit{max n} (high) rather than depicting a probability in this column. (\textbf{right}) The median migratory probability of experiencing above threshold temperatures.}
\label{fig:5}
\end{figure}

Reflecting regionally elevated thermal exposure probabilities, cummulative migratory thermal exposure generally increased with migratory distance and was much higher in the South Thompson River watershed (Figs. $\ref{fig:5}$, $\ref{fig:6}$, S3). In other words, the opportunity to experience highly likely (>50\%) stressful temperatures was greater for those migrating further along the mainstem of the South Thompson River watershed. This finding was particularly true in 2014, early and mid-summer of 2015 and mid-summer of 2016, while thermal exposure was relatively unlikely in 2017. Cummulative migratory thermal exposure above 22$\text{\textdegree}$C, was nearly non-existent in either basin across sites, years and migratory periods (Fig. S4). 

Contrary to the cummulative thermal exposure, the cummulative probability of thermal exposure declined with migratory distance in the South Thompson (Figs. $\ref{fig:5}$, S5). In the North Thompson, the probability of thermal stress was notably highest in the western portion of the watershed and did not follow the migratory distance patterns of the South Thompson. These findings also reflect the thermal exposure probability dynamics of the mainstem versus the tributaries (Fig. $\ref{fig:5}$). Moving further up the watershed and away from the mainstem reduces the local probability of thermal stress resulting in declining overall migratory probabilities of thermal stress. There was little to no evidence that the 22$\text{\textdegree}$C threshold is likely among any migration route (Fig. S6).

*Chinook Salmon Thermal Exsposure*

The average thermal exposure for Thompson River watershed chinook populations demonstrates that South Thompson River watershed populations are much more likely to experience thermal stress than those in the North Thompson River watershed (Fig. $\ref{fig:6}$). In many years the potential for thermal stress (>19$\text{\textdegree}$C) among southern populations persisted after diverging off the mainstem and arriving on the spawning grounds but never above 50\% probability. At the 22$\text{\textdegree}$C threshold, only one population exhibited any likely (>50\%) chance of experiencing these temperatures during migration and they were extremely unlikely across populations on the spawning grounds (Fig. $\ref{fig:6}$). Notably, the number of likely stressful temperatures declined with shorter migration routes but spawn site thermal exposure potential increased. Across populations, thermal exposure potential declined from 2014 to 2017.

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{./Water_Temperature/images/chinook.png}
	\caption{Thermal exposure above 19 and 22\text{\textdegree}C for 16 chinook salmon populations in the North (circles) and South (triangles) Thompson River watersheds from 2014 to 2017. (\textbf{left}) The cummulative number of likely (>50\%) above threshold mean temperature estimates during migration. (\textbf{right}) The mean probability of experiencing above threshold temperatures upon arriving on the spawning grounds. All values are calculated given the populations expected migratory period and populations are ordered from furthest to shortest migration within their watershed.}
\label{fig:6}
\end{figure}

#Discussion

Persistent across measures of migratory thermal exposure is the tendency for lake dominated regions to result in higher levels of stress potential. The positive effect of lakes, quantified in the *A* and $\alpha$ coefficient SNNM models, are clearly observed in the South Thompson watershed and western portion of the North Thompson watershed. Interestingly, glacial and high elevation lake effects dampened thermal stress potential but these landscape variables are limited on the landscape and may explain lower thermal exposure potential in the North Thompson River watershed relative to the South Thompson where these landscape attributes are non-exisitent.

Southern system seems like a small increase in temperature is likely to persistently lead to stressful thermal conditions whereas the north Thompson has large regions buffered from the warming climate.

El Nino and PDO leading to higher levels of winter snowfall must be playing a significant role but is currently not being considered. This may explain the unexplained variance in the hysteresis model.

Meta-population stability? Identifying at risk population because typically thermal stress is only considered on the spawning ground rather than the entire migratory route.

The model does not capture variation around the mean very well but it's important to note that these are point source estimates of temperature and may not reflect the greater heterogeneity of the stream that can be exploited by a salmon. This is why we need to consider what is possible considering 95% of variance around the mean.

I've noticed that at times the temporal model underestimates parameters (e.g., alpha) but the spatial model often is able to correct for these underestimates. It would seem this is because the landscape relationships result in parameter expectations built on largely accurate coefficient estimates that correct for the occasional underestimate. On the other hand, there is an underestimate that is further under-estimated by the spatial model on the Adams river. This may be because ...

\newpage

#References
