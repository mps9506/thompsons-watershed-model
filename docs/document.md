---
title: "Empirical Streamflow Predictions at Thompsons Creek"
date: "2021-03-31"
github-repo: https://github.com/mps9506/thompson-stage-discharge
bibliography: bibliography.bib
biblio-style: "apalike"
link-citations: true
---




```r
##note, need development version of bookdown and gt to properly compile this document

## readr imports data
library(readr)
library(tidyr)
## tibbles are advanced fancy dataframes
library(tibble)
## dplyr for data handling and piping functions
library(dplyr)
## ggplot for plots
library(ggplot2)
## stringr to read and manipulate strings
library(stringr)
## here is a function to ensure file paths are correct
library(here)
## units and ggforce facilitate attaching units to data
library(units)
library(ggforce)
## hrbrtheme is optional, I use it to pretty my plots
library(hrbrthemes)
## patchwork and cowplot support arranging multiple ggplots into one plot
library(patchwork)
library(cowplot)
## lubridate provides functions for handling time and date
library(lubridate)
## purrr lets us use map_ functions as an alternative to loops
library(purrr)
## hydroGOF provide goodness of fit metrics (NSE, RMSE, etc.)
library(hydroGOF)
## tsibble and imputeTS will allow some simple time series interpolation
library(tsibble)
library(imputeTS)
## gt
library(gtsummary)
library(flextable)
## nls.multstart fits non-linear least squares using the Levenberg-Marquardt algorithm with multiple starting values.
library(nls.multstart)
##
library(mgcv)
library(gratia)
## apply drainage area ratio
library(dartx)

## set some options
update_geom_font_defaults(font_rc)
units_options(parse = FALSE)


## some custom functions
theme_ms <- function(...) {
  theme_ipsum_rc(plot_margin = margin(10,10,10,10),
              axis_title_just = "c") +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = "transparent", 
            colour = NA), 
          panel.border = element_rect(fill = NA, 
            colour = "grey20"),
          plot.title = element_text(size = 11,
                                    face = "plain"),
          ...)
}

exponent <- function(x, pow) {
  (abs(x)^pow)*sign(x)
}
```

## About

This is an exploratory document investigating the feasibility of estimating mean daily streamflows in the Thompsons Creek watershed using empirical methods. In a previous document mean daily streamflows at three sites in the watershed were developed for March 2020 through March 20221 using measured flows, depths, and rating curves. The one-year record of streamflows will be used to fit and validate drainage area ratio, regression, and semi-parametric regression based methods for estimating streamflow using locally available data.

## Introduction


 


### Ungaged streamflow estimation

**Statistical information transfer methods**

Statistical information transfer and empirical regression are two relatively simple methods for estimating streamflows in poorly gaged watersheds using only measured stremaflow values from nearby gaged watersheds. Statistical transfer procedures simply transfer flow duration curves or daily streamflow values from a gaged watershed to the ungaged watershed using assumed relationships between area and runoff. The most common statistical transfer method is the drainage area ratio. With the drainage area ratio, daily streamflows are transferred from one basin to the other by multiplying the area ratio to daily streamflows:

\begin{equation}
Q_y^t = Q_x^t\bigg(\frac{A_y}{A_x}\bigg)^\phi
  (\#eq:dar)
\end{equation}

Where $Q_y^t$ is streamflow at ungaged basin $y$ and time $t$, $Q_x^t$ is streamflow at gaged basin $x$ and time $t$, and $\frac{A_y}{A_x}$ is the area ratio of the basins. Parameter $\phi$ is typically equal to one  [@asquith_statewide_2006]. However, @asquith_statewide_2006 provides empirically estimated values of $\phi$ for use in the drainage area ratio when applied in Texas. With an available short term streamflow record available, various stream gages can be assessed for performance using the drainage area ratio method. However, when that short-term streamflow record is available we extend the simple drainage area ratio to a linear regression for streamflow estimation using one or more gaged watersheds:

\begin{equation}
Q_y = \beta_0 + \beta_n{Q_{xn}} + \varepsilon
  (\#eq:linearregression1)
\end{equation}

where $Q_y$ is the predicted mean daily streamflow at the ungaged or temporarily gaged site, $\beta_0$ is the intercept, $Q_{xn}$ are mean daily streamflows at gaged watershed $n$, $\beta_n$ is the regression coefficient, and $varepsilon$ is the residual error term assumed normally distributed around mean zero. A linear regression of this form still acts as a streamflow transfer methods like the drainage area approach but allows for an easy incorporation of additional model terms such as lagged streamflows which might improve predictive performance.

**Semi-parametric rainfall-runoff regression**

If nearby gaged watersheds are unavailable or are not reflective of the streamflow responses in the ungaged watershed, locally available weather information can be used to empircally estimate streamflow response. A number of empirically based rainfall-runoff routing models are available that account for soil and land-use conditions to predict streamflow response (SIMHYD, IHACRES, and Sacramento rainfall-runoff models are examples). More complex mechanistic models that simulate hydrologic and water quality responses to landuse and precipitation are also available (SWAT and HSPF are two examples). The mechanistic models have a steep requirement for data and ability of the technician developing the model. The routing models and mechanistic models are outside the scope of work for this particular project. Here we focus on using a semi-parametric regression based approach to predict streamflow using locally available weather data. Empirical regression based approaches require a period of measured streamflow and some predictor variables to estimate the runoff response from. Typically, daily rainfall and temperature data are employed to fit a regression model to measured streamflow response:

\begin{equation}
Q_i = \beta_0 + \beta_1x_i+ \varepsilon_i
  (\#eq:linearregression2)
\end{equation}

where $Q$ is predicted discharge on day $i$, $\beta_0$ is the intercept, $\beta_1$ is the regression coefficient, $x$ is the predictor variable value on the $i$th day. The error term, $\varepsilon_i$ is assumed normally distributed around mean zero. Using simple linear or multiple linear regression $Q$ is typically log transformed prior to estimating the regression coefficients. Generalized linear models (GLMs) can instead be used when the data distribution is skewed (more specifically, when the residuals are not expected to be normally distributed). This is accomplished through the inclusion of a link function that describes how the mean of the response depends on the linear predictor and a variance function the describes how variance depends on the mean. Using the GLM and appropriate link function, transformation of the response variable can be avoided. This approach is highly desirable due to biases in the variance structure that are introduced when estimating the means of log transformed data that is intended to be back transformed to the original scale.

If the relationship between predictor variables and the response are expected to be nonlinear, polynomial terms can be included. However, an extension of GLMs called generalized additive models (GAMs) allows relatively easy fitting of these nonlinear terms to the data. With GAMs, the response variable depends on the sum of smoothing functions applied to each predictor variable:

\begin{equation}
Q_i = \beta_0 + f(x_1) + \varepsilon
  (\#eq:gam)
\end{equation}

where $Q$ is predicted discharge on day $i$, $\beta_0$ is the intercept, and some function $f(x_1)$ is the linear predictor. In the case of GAMs fit using the `mgcv` package in R, $f$ is a smoothing function fit to the data using generalized cross validation or restricted maximum likelihood. Due to the smoothing function and link functions, GAMs are extremely flexible for fitting regression models to data of different distributions and responses. However, compared to linear regression and GLMs the inclusion of the smoothing functions limits interpretability because traditional regression coefficients are not part of the model structure. Therefore, the effect of each individual smoothing function on the response variable mean is shown graphically.



## Method

### Data

Two National Oceanic and Atmospheric Administration (NOAA) Global Historical Climatology Network (GHCN) locations provide daily precipitation data for the project area (Figure \@ref(fig:precipsum); Figure \@ref(fig:temp)). GHCND daily summaries were downloaded for GHCND:USW00003904 (Easterwood Airport) using the NOAA API services and the `rnoaa` package in R [@chamberlain_2019].


```r
easterwood_precip <- read_csv("Data/noaa_precip/easterwood_precip.csv",
                               col_types = cols(
                                 date = col_datetime(format = ""),
                                 datatype = col_character(),
                                 station = col_character(),
                                 value = col_double(),
                                 fl_m = col_character(),
                                 fl_q = col_character(),
                                 fl_so = col_character(),
                                 fl_t = col_character(),
                                 units = col_character())) %>%
  mutate(value =  set_units(value/10, "mm")) %>%
  select(date, datatype, station, value)
  


easterwood_precip %>%
  mutate(date = as.Date(date),
         value = as.numeric(value)) %>%
  ggplot() +
  geom_line(aes(date, value)) +
  geom_point(aes(date, value), alpha = 0) +
  labs(y = "Daily precipitation (mm)",
       x = "Date") +
  facet_wrap(~station, ncol = 1) +
  theme_ms(plot.margin = margin(10,0,10,10),
           panel.spacing = unit(0, "lines"),
           panel.spacing.x = unit(0, "lines"),
           panel.spacing.y = unit(2, "lines")) -> p1

  

easterwood_precip %>%
  mutate(date = as.Date(date),
         value = as.numeric(value)) %>%
  ggplot() +
  geom_histogram(aes(value), binwidth = 5) +
  facet_wrap(~station, ncol = 1) +
  coord_flip() +
  theme_ms() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(10,10,10,0),
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(2, "lines"),
        strip.text = element_text(color="transparent")) -> p2


plot_grid(p1,p2,align = "h", axis = "bt", rel_widths = c(2, 1))
```

<div class="figure">
<img src="document_files/figure-html/precipsum-1.png" alt="Daily precipitation and marginal histograms of daily precipitation." width="100%" />
<p class="caption">(\#fig:precipsum)Daily precipitation and marginal histograms of daily precipitation.</p>
</div>



```r
easterwood_tmax <- read_csv("Data/noaa_precip/easterwood_tmax.csv",
                               col_types = cols(
                                 date = col_datetime(format = ""),
                                 datatype = col_character(),
                                 station = col_character(),
                                 value = col_double(),
                                 fl_m = col_character(),
                                 fl_q = col_character(),
                                 fl_so = col_character(),
                                 fl_t = col_character(),
                                 units = col_character())) %>%
  mutate(value =  set_units(value/10, "°C")) %>%
  select(date, datatype, station, value)

easterwood_tmax %>%
  mutate(date = as.Date(date),
         value = as.numeric(value)) %>%
  ggplot() +
  geom_step(aes(date, value), alpha = 0.5) +
  labs(y = "Daily Maximum Temperature [°C]",
       x = "Date") +
  theme_ms(plot.margin = margin(10,0,10,10),
           panel.spacing = unit(0, "lines"),
           panel.spacing.x = unit(0, "lines"),
           panel.spacing.y = unit(2, "lines")) -> p1

easterwood_tmax %>%
  mutate(date = as.Date(date),
         value = as.numeric(value)) %>%
  ggplot() +
  geom_density(aes(value), 
                 binwidth = 1.1,
                 alpha = 0.5) +
  coord_flip() +
  theme_ms() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(10,10,10,0),
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(2, "lines"),
        strip.text = element_text(color="transparent")) -> p2

plot_grid(p1,p2,align = "h", axis = "bt", rel_widths = c(2, 1))
```

<div class="figure">
<img src="document_files/figure-html/temp-1.png" alt="Daily maximum temperature and marginal histograms of daily maximum temperature." width="100%" />
<p class="caption">(\#fig:temp)Daily maximum temperature and marginal histograms of daily maximum temperature.</p>
</div>

We identified two wastewater treatment facilities located upstream of the lowest discharge point. The Still Creek WWTF (TPDES Permit No. WQ0010426002) is permitted to discharge 4.0 MGD to Still Creek and discharge flows past SWQM sites 16882 and 16396. The Sanderson Farm Inc. facility (TPDES Permit No. WQ0003821000) is permitted to discharge 1.678 MGD to a tributary of Cottonwood Branch and the discharge flows past SWQM site 16396. Mean daily wastewater facility discharges were downloaded from the EPA ECHO database (Figure \@ref(fig:wwtf)) using the `echor` R package and subtracted from the measured mean daily streamflows to better represent naturalized mean daily flows [@schramm_2020].


```r
wwtf <- read_csv("Data/EPA_WWTF/mean_daily_discharges.csv",
                 col_types = cols(
                   npdes_id = col_character(),
                   date = col_date(format = ""),
                   mgd = col_double(),
                   cfs = col_double()))
wwtf %>%
  ggplot() +
  geom_step(aes(date, cfs, color = npdes_id)) +
  scale_x_date(date_breaks = "month",
               date_labels = "%b-%y") +
  scale_color_discrete(labels = c("WQ0010426002, Still Creek WWTF",
                                "WQ0003821000, Sanderson Farm Inc." )) +
  labs(x = "Date", y = "Mean Daily Discharge [cfs]") + 
  theme_ms() +
  guides(x = guide_axis(angle = 90)) +
  theme(axis.line.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.length = grid::unit(5, "pt"),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())
```

<div class="figure">
<img src="document_files/figure-html/wwtf-1.png" alt="Mean daily WWTF discharges" width="100%" />
<p class="caption">(\#fig:wwtf)Mean daily WWTF discharges</p>
</div>


### Daily flow estimation

Mean daily streamflows from 2020-03-03 through 2021-01-14 were estimated using rating curves at SWQM stations 16396, 16397, and 16882 (Table \@ref(tab:dartable)). Mean daily wastewater discharges were subtracted from the flow record to represent naturalized flows. Mean daily temperature and total daily precipitation were joined to the data record by date. Figure \@ref(fig:meandailyresults) shows the hydrograph of naturalized flows and precipitation at each SWQM station.

Daily naturalized flow was estimated using the drainage area ratio equation (Equation \@ref(eq:dar)) at USGS stream gages 08065800, 08109800, and 08110100 (Table \@ref(tab:dartable)). Values of \phi recommended in @asquith_statewide_2006 were used. Linear regression was also used to estimate log transformed flows at each SWQM station using mean daily flows and 1-day lagged mean daily flows from each USGS stream gage as a predictor variable.


```r
dar_table <- tibble(Site = c("SWQM-16396", "SWQM-16397", "SWQM-16882", "USGS-08065800", "USGS-08109800", "USGS-08110100"),
       Description = c("Thompsons Creek at Silver Hill Rd",
                       "Thompsons Creek at Hwy 21",
                       "Still Creek at Hwy 21",
                       "Bedias Creek near Madisonville",
                       "East Yegua Creek near Dime Box",
                       "Davidson Creek near Lyons"),
       Area = c(42.33,24.21,10.03,321,244,195))

kable(dar_table,
      caption = "TCEQ SWQM stations and USGS stream gages used to develop flows with drainage area ratio and linear regression methods.")
```



Table: (\#tab:dartable)TCEQ SWQM stations and USGS stream gages used to develop flows with drainage area ratio and linear regression methods.

|Site          |Description                       |   Area|
|:-------------|:---------------------------------|------:|
|SWQM-16396    |Thompsons Creek at Silver Hill Rd |  42.33|
|SWQM-16397    |Thompsons Creek at Hwy 21         |  24.21|
|SWQM-16882    |Still Creek at Hwy 21             |  10.03|
|USGS-08065800 |Bedias Creek near Madisonville    | 321.00|
|USGS-08109800 |East Yegua Creek near Dime Box    | 244.00|
|USGS-08110100 |Davidson Creek near Lyons         | 195.00|





```r
##naturalize streamflows

df <- read_csv("Data/daily_streamflow/model_df.csv",
               col_types = cols(
                 Site = col_character(),
                 date = col_date(format = ""),
                 mean_daily = col_double()
                 ))


df %>%
  bind_rows(tibble(Site = c(rep("16396",5),rep("16397",5),rep("16882",5)),
                   date = rep(seq.Date(as.Date("2020-02-27"), as.Date("2020-03-02"), by = "day"),3),
                   mean_daily = NA)) %>%
  arrange(date) %>%
  left_join(easterwood_precip, by = c("date" = "date")) %>%
  mutate(value = as.numeric(value)) %>%
  dplyr::rename(ewood_precip = value) %>%
  left_join(easterwood_tmax, by = c("date" = "date")) %>%
  mutate(value = as.numeric(value)) %>%
  dplyr::rename(ewood_tmax = value) %>%
  dplyr::select(Site, date, mean_daily, ewood_precip, ewood_tmax) %>%
  left_join(wwtf %>% pivot_wider(id_cols = date, names_from = npdes_id, values_from = cfs),
            by = c("date" = "date")) %>%
  ## remove WWTF influence from discharge record
  mutate(adjusted_flow = case_when(
    Site == 16396 ~ mean_daily - TX0025071 - TX0113603, 
    Site == 16882 ~ mean_daily - TX0025071,
    Site == 16397 ~ mean_daily)) %>%
  mutate(adjusted_flow = case_when(
    adjusted_flow < 0 ~ 0,
    adjusted_flow >= 0 ~ adjusted_flow)) %>%
  mutate(doy = lubridate::yday(date))-> df
```



```r
## plot hydrograph of Thompsons @ Silver Hill
p1 <- ggplot(df %>% filter(Site == "16396") %>% mutate(date = as.Date(date))) + 
  geom_line(aes(date, adjusted_flow, color = "Mean Daily Streamflow")) +
  scale_y_continuous(position = "left", 
                  limits = c(0,600),
                  expand = c(0,0)) +
  scale_x_date(date_breaks = "month",
               date_labels = "%b-%y") +
  labs(y = "Mean daily streamflow [cfs]", 
       x = "Date",
       caption = "16396, Thompsons Creek at Silver Hill Rd") +
  scale_color_manual(values = c("dodgerblue")) +
  theme_ms() +
  guides(x = guide_axis(angle = 90)) +
  theme(axis.text.x = element_text(size = 8),
        axis.title.y.left = element_text(hjust = 0),
        axis.ticks.y.left = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.length = grid::unit(5, "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none") 

p2 <- ggplot(df %>% filter(Site == "16396")) + 
  geom_line(aes(date, ewood_precip, color = "Total Daily Preciptitation")) +
  scale_y_reverse(position = "right", 
                  limits = c(350,0),
                  breaks = c(0,25,50),
                  labels = c(0,25,50),
                  expand = c(0,0)) +
  labs(y = "Total Daily Precipitaiton [mm]") +
  scale_color_manual(values = c("dodgerblue4")) +
  theme_ms() +
  guides(x = guide_axis(angle = 90)) +
  theme(axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y.right = element_text(hjust = 1),
    axis.title.x = element_blank(),
    axis.ticks.y.right = element_line(color = "black"),
    axis.ticks.length = grid::unit(5, "pt"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
    ) 
set_null_device("agg")
aligned_plots <- align_plots(p1, p2, align="hv", axis="tblr")
set_null_device("agg")
hg_plot1 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

## plot hydrograph of Thompsons @ Hwy21
p1 <- ggplot(df %>% filter(Site == "16397") %>% mutate(date = as.Date(date))) + 
  geom_line(aes(date, adjusted_flow, color = "Mean Daily Streamflow")) +
  scale_y_continuous(position = "left", 
                  limits = c(0,75),
                  expand = c(0,0)) +
  scale_x_date(date_breaks = "month",
               date_labels = "%b-%y") +
  labs(y = "Mean daily streamflow [cfs]", 
       x = "Date",
       caption = "16397, Thompsons Creek at Hwy 21") +
  scale_color_manual(values = c("dodgerblue")) +
  theme_ms() +
  guides(x = guide_axis(angle = 90)) +
  theme(axis.text.x = element_text(size = 8),
        axis.title.y.left = element_text(hjust = 0),
        axis.ticks.y.left = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.length = grid::unit(5, "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none") 

p2 <- ggplot(df %>% filter(Site == "16397")) + 
  geom_line(aes(date, ewood_precip, color = "Total Daily Preciptitation")) +
  scale_y_reverse(position = "right", 
                  limits = c(350,0),
                  breaks = c(0,25,50),
                  labels = c(0,25,50),
                  expand = c(0,0)) +
  labs(y = "Total Daily Precipitaiton [mm]") +
  scale_color_manual(values = c("dodgerblue4")) +
  theme_ms() +
  guides(x = guide_axis(angle = 90)) +
  theme(axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y.right = element_text(hjust = 1),
    axis.title.x = element_blank(),
    axis.ticks.y.right = element_line(color = "black"),
    axis.ticks.length = grid::unit(5, "pt"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
    ) 
set_null_device("agg")
aligned_plots <- align_plots(p1, p2, align="hv", axis="tblr")
set_null_device("agg")
hg_plot2 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])


## plot hydrograph of Sill Creek @ Hwy21
p1 <- ggplot(df %>% filter(Site == "16882") %>% mutate(date = as.Date(date))) + 
  geom_line(aes(date, adjusted_flow, color = "Mean Daily Streamflow")) +
  scale_y_continuous(position = "left", 
                  limits = c(0,75),
                  expand = c(0,0)) +
  scale_x_date(date_breaks = "month",
               date_labels = "%b-%y") +
  labs(y = "Mean daily streamflow [cfs]", 
       x = "Date",
       caption = "16882, Sill Creek at Hwy 21") +
  scale_color_manual(values = c("dodgerblue")) +
  theme_ms() +
  guides(x = guide_axis(angle = 90)) +
  theme(axis.text.x = element_text(size = 8),
        axis.title.y.left = element_text(hjust = 0),
        axis.ticks.y.left = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.length = grid::unit(5, "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none") 

p2 <- ggplot(df %>% filter(Site == "16882")) + 
  geom_line(aes(date, ewood_precip, color = "Total Daily Preciptitation")) +
  scale_y_reverse(position = "right", 
                  limits = c(350,0),
                  breaks = c(0,25,50),
                  labels = c(0,25,50),
                  expand = c(0,0)) +
  labs(y = "Total Daily Precipitaiton [mm]") +
  scale_color_manual(values = c("dodgerblue4")) +
  theme_ms() +
  guides(x = guide_axis(angle = 90)) +
  theme(axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y.right = element_text(hjust = 1),
    axis.title.x = element_blank(),
    axis.ticks.y.right = element_line(color = "black"),
    axis.ticks.length = grid::unit(5, "pt"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
    ) 
set_null_device("agg")
aligned_plots <- align_plots(p1, p2, align="hv", axis="tblr")
set_null_device("agg")
hg_plot3 <- ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])


hg_plot1 / hg_plot2 / hg_plot3
```

<div class="figure">
<img src="document_files/figure-html/meandailyresults-1.png" alt="Hydrograph of mean daily streamflows at SWQM 16396, 16397, and 16882 with reported total daily precipitation at Easterwood Airport." width="100%" />
<p class="caption">(\#fig:meandailyresults)Hydrograph of mean daily streamflows at SWQM 16396, 16397, and 16882 with reported total daily precipitation at Easterwood Airport.</p>
</div>

GLMs and GAMs relied on local precipitation and temperature data to estimate streamflows at each SWQM station (Figure \@ref(fig:precipsum); Figure \@ref(fig:temp)). The GLM was of form:

\begin{multline}
Q = P + T + DOY + P_{sum,3} + T_{mean,5} + (P\times T) + (P \times P_{sum,3})\\ + (P \times T_{mean,5}) + (T \times P_{sum,3})\\ + (T \times T_{mean,5}) + (P_{sum,3} \times T_{mean,5})
  (\#eq:glmstructure)
\end{multline}

where $P$ is total daily precipitation and $T$ is mean daily temperature, and are assumed to be the main forcing variables influencing streamflow. $DOY$ is day of the year and included as a seasonal predictor. $P_{sum,3}$ is the rolling 3-day sum rainfall and included as an indicator of wetness in the watershed. $T_{mean,5}$ is the rolling 5-day mean temperature and included as an indicator of seasonal temperature condition that is less daily variance than $T$ and maybe an improved indicator of potential evapotransportation conditions. The GLM error structure was fit with a scaled t distribution and a log link function.

The GAM included the same variables but with smoothing functions:

\begin{multline}
Q = f(P) + f(T) + f(DOY) + f(P_{sum,3}) + f(T_{mean,5})\\ + f(P,T) + f(P, P_{sum,3}) + f(P,T_{mean,5})\\ + f(T,P_{sum,3}) + (T,T_{mean,5}) + f(P_{sum,3},T_{mean,5})
  (\#eq:gamstructure)
\end{multline}

where $f()$ is the smoothing function. Where a single covariate is smoothed, a thin plate regression spline was fit to the data using restricted maximum likelihood to estimated the automatically select the optimal smoothing parameters. For smoothing functions applied to two parameters, a tensor product smooth function was applied. The GAM error structure was fit using the scaled t distribution and a log link function.

Performance of predicted streamflows over the period of record for each method was assessed using Nash-Sutcliffe Efficiency (NSE) and Kling-Gupta Efficiency (KGE) goodness-of-fit metrics. 

## Results





```r
## create a model df for 16396
## additional variables = 
## lagPrecip = 1 day lag precipitation
## wetness index = 3 day total precip
## ET index = 5 day avg tmax


df %>%
  dplyr::filter(Site == "16396") %>%
  arrange(date) -> df_16396

df_16396 %>%
  mutate(lagPrecip = lag(ewood_precip),
         wetness = map(row_number(.$date),
                    ~{if(.x - 1 <= 0) {df_16396$ewood_precip[.x]}
                      if(.x - 2 <= 0) {sum(df_16396$ewood_precip[.x],
                                           df_16396$ewood_precip[.x-1])}
                      if(.x - 2 > 0) {
                       sum(
                        df_16396$ewood_precip[.x],
                        df_16396$ewood_precip[.x-1],
                        df_16396$ewood_precip[.x-2],
                        na.rm = TRUE
                        ) 
                      }}),
         et = map(row_number(.$date),
                    ~{if(.x - 1 <= 0) {df_16396$ewood_tmax[.x]}
                      if(.x - 2 <= 0) {df_16396$ewood_tmax[.x-1]}
                      if(.x - 3 <= 0) { mean(c(df_16396$ewood_tmax[.x-1],
                                            df_16396$ewood_tmax[.x-2]),
                                            na.rm = TRUE)}
                      if(.x - 4 <= 0) { mean(c(df_16396$ewood_tmax[.x-1],
                                            df_16396$ewood_tmax[.x-2],
                                            df_16396$ewood_tmax[.x-3]),
                                            na.rm = TRUE)}
                      if(.x - 5 <= 0) { mean(c(df_16396$ewood_tmax[.x-1],
                                            df_16396$ewood_tmax[.x-2],
                                            df_16396$ewood_tmax[.x-3],
                                            df_16396$ewood_tmax[.x-4]),
                                            na.rm = TRUE)}
                      mean(c(df_16396$ewood_tmax[.x-1],
                        df_16396$ewood_tmax[.x-2],
                        df_16396$ewood_tmax[.x-3],
                        df_16396$ewood_tmax[.x-4],
                        df_16396$ewood_tmax[.x-5]),
                        na.rm = TRUE)})) %>%
  unnest(c(wetness, et)) %>%
  dplyr::filter(!is.na(mean_daily)) -> df_16396


## create a model df for 16397
## additional variable = 
## wetness index = 3 day total precip
## ET index = 5 day avg tmax

df %>%
  dplyr::filter(Site == "16397") %>%
  arrange(date) -> df_16397

df_16397 %>%
  mutate(lagPrecip = lag(ewood_precip),
         wetness = map(row_number(.$date),
                    ~{if(.x - 1 <= 0) {df_16397$ewood_precip[.x]}
                      if(.x - 2 <= 0) {sum(df_16397$ewood_precip[.x],
                                           df_16397$ewood_precip[.x-1])}
                      if(.x - 2 > 0) {
                       sum(
                        df_16397$ewood_precip[.x],
                        df_16397$ewood_precip[.x-1],
                        df_16397$ewood_precip[.x-2],
                        na.rm = TRUE
                        ) 
                      }}),
         et = map(row_number(.$date),
                    ~{if(.x - 1 <= 0) {df_16397$ewood_tmax[.x]}
                      if(.x - 2 <= 0) {df_16397$ewood_tmax[.x-1]}
                      if(.x - 3 <= 0) { mean(c(df_16397$ewood_tmax[.x-1],
                                            df_16397$ewood_tmax[.x-2]),
                                            na.rm = TRUE)}
                      if(.x - 4 <= 0) { mean(c(df_16397$ewood_tmax[.x-1],
                                            df_16397$ewood_tmax[.x-2],
                                            df_16397$ewood_tmax[.x-3]),
                                            na.rm = TRUE)}
                      if(.x - 5 <= 0) { mean(c(df_16397$ewood_tmax[.x-1],
                                            df_16397$ewood_tmax[.x-2],
                                            df_16397$ewood_tmax[.x-3],
                                            df_16397$ewood_tmax[.x-4]),
                                            na.rm = TRUE)}
                      mean(c(df_16397$ewood_tmax[.x-1],
                        df_16397$ewood_tmax[.x-2],
                        df_16397$ewood_tmax[.x-3],
                        df_16397$ewood_tmax[.x-4],
                        df_16397$ewood_tmax[.x-5]),
                        na.rm = TRUE)})) %>%
  unnest(c(wetness, et)) %>%
  dplyr::filter(!is.na(mean_daily)) -> df_16397

## create a model df for 16882
## additional variable = 
## wetness index = 3 day total precip
## ET index = 5 day avg tmax

df %>%
  dplyr::filter(Site == "16882") %>%
  arrange(date) -> df_16882

df_16882 %>%
  mutate(lagPrecip = lag(ewood_precip),
         wetness = map(row_number(.$date),
                    ~{if(.x - 1 <= 0) {df_16882$ewood_precip[.x]}
                      if(.x - 2 <= 0) {sum(df_16882$ewood_precip[.x],
                                           df_16882$ewood_precip[.x-1])}
                      if(.x - 2 > 0) {
                       sum(
                        df_16882$ewood_precip[.x],
                        df_16882$ewood_precip[.x-1],
                        df_16882$ewood_precip[.x-2],
                        na.rm = TRUE
                        ) 
                      }}),
         et = map(row_number(.$date),
                    ~{if(.x - 1 <= 0) {df_16882$ewood_tmax[.x]}
                      if(.x - 2 <= 0) {df_16882$ewood_tmax[.x-1]}
                      if(.x - 3 <= 0) { mean(c(df_16882$ewood_tmax[.x-1],
                                            df_16882$ewood_tmax[.x-2]),
                                            na.rm = TRUE)}
                      if(.x - 4 <= 0) { mean(c(df_16882$ewood_tmax[.x-1],
                                            df_16882$ewood_tmax[.x-2],
                                            df_16882$ewood_tmax[.x-3]),
                                            na.rm = TRUE)}
                      if(.x - 5 <= 0) { mean(c(df_16882$ewood_tmax[.x-1],
                                            df_16882$ewood_tmax[.x-2],
                                            df_16882$ewood_tmax[.x-3],
                                            df_16882$ewood_tmax[.x-4]),
                                            na.rm = TRUE)}
                      mean(c(df_16882$ewood_tmax[.x-1],
                        df_16882$ewood_tmax[.x-2],
                        df_16882$ewood_tmax[.x-3],
                        df_16882$ewood_tmax[.x-4],
                        df_16882$ewood_tmax[.x-5]),
                        na.rm = TRUE)})) %>%
  unnest(c(wetness, et)) %>%
  dplyr::filter(!is.na(mean_daily)) %>%
  mutate(non_zero = case_when(
    adjusted_flow == 0 ~ 0,
    adjusted_flow > 0 ~ 1
  ))-> df_16882
```

### Statistical Information Transfer

**DAR**

Mean daily streamflows estimated using DAR applied to three different USGS streamgages are displayed in Figure \@ref[fig:dar16396results]. All three USGS gages results in biased results at low streamflows (consistently underpredicted streamflows). The hydrographs indicate the timing and magnitude of stormflow events are routinely mismatched using any of the three USGS gages. This visual validation suggests DAR method is not suitable in the Thompsons Creek watershed.


```r
## apply DAR to 16396 using each select USGS gage
usgs_08065800 <- readr::read_csv("Data/USGS_Streamflow/08065800.csv") %>%
  dplyr::rename(Flow_08065800 = "133946_00060_00003") %>%
  dplyr::select(datetime, Flow_08065800)

usgs_08109800 <- readr::read_csv("Data/USGS_Streamflow/08109800.csv") %>%
  dplyr::rename(Flow_08109800 = `135356_00060_00003`) %>%
  dplyr::select(datetime, Flow_08109800)

usgs_08110100 <- readr::read_csv("Data/USGS_Streamflow/08110100.csv") %>%
  dplyr::rename(Flow_08110100 = `135389_00060_00003`) %>%
  dplyr::select(datetime, Flow_08110100)


df_16396 %>%
  dplyr::select(Site, date, adjusted_flow) %>%
  left_join(usgs_08065800, by = c("date" = "datetime")) %>%
  dartx(Flow_08065800, 42.33/321) %>%
  dplyr::rename(DAR_Q_08065800 = Q,
                Q_percentile_08065800 = Q_percentile,
                exp_08065800 = exp) %>%
  left_join(usgs_08109800,  by = c("date" = "datetime")) %>%
  dartx(Flow_08109800, 42.33/244) %>%
  dplyr::rename(DAR_Q_08109800 = Q,
                Q_percentile_08109800 = Q_percentile,
                exp_08109800 = exp) %>%
  left_join(usgs_08110100,  by = c("date" = "datetime")) %>%
  dartx(Flow_08110100, 42.33/195) %>%
  dplyr::rename(DAR_Q_08110100 = Q,
                Q_percentile_08110100 = Q_percentile,
                exp_08110100 = exp) -> dar_results_16396

df_16882 %>%
  dplyr::select(Site, date, adjusted_flow) %>%
  left_join(usgs_08065800, by = c("date" = "datetime")) %>%
  dartx(Flow_08065800, 24.21/321) %>%
  dplyr::rename(DAR_Q_08065800 = Q,
                Q_percentile_08065800 = Q_percentile,
                exp_08065800 = exp) %>%
  left_join(usgs_08109800,  by = c("date" = "datetime")) %>%
  dartx(Flow_08109800, 24.21/244) %>%
  dplyr::rename(DAR_Q_08109800 = Q,
                Q_percentile_08109800 = Q_percentile,
                exp_08109800 = exp) %>%
  left_join(usgs_08110100,  by = c("date" = "datetime")) %>%
  dartx(Flow_08110100, 24.21/195) %>%
  dplyr::rename(DAR_Q_08110100 = Q,
                Q_percentile_08110100 = Q_percentile,
                exp_08110100 = exp) -> dar_results_16882

df_16397 %>%
  dplyr::select(Site, date, adjusted_flow) %>%
  left_join(usgs_08065800, by = c("date" = "datetime")) %>%
  dartx(Flow_08065800, 10.03/321) %>%
  dplyr::rename(DAR_Q_08065800 = Q,
                Q_percentile_08065800 = Q_percentile,
                exp_08065800 = exp) %>%
  left_join(usgs_08109800,  by = c("date" = "datetime")) %>%
  dartx(Flow_08109800, 10.03/244) %>%
  dplyr::rename(DAR_Q_08109800 = Q,
                Q_percentile_08109800 = Q_percentile,
                exp_08109800 = exp) %>%
  left_join(usgs_08110100,  by = c("date" = "datetime")) %>%
  dartx(Flow_08110100, 10.03/195) %>%
  dplyr::rename(DAR_Q_08110100 = Q,
                Q_percentile_08110100 = Q_percentile,
                exp_08110100 = exp) -> dar_results_16397
```


```r
dar_results_16396 %>%
  dplyr::select(Site, date, adjusted_flow, DAR_Q_08065800, DAR_Q_08109800, DAR_Q_08110100) %>%
  pivot_longer(c(DAR_Q_08065800, DAR_Q_08109800, DAR_Q_08110100),
               names_to = "Source_Site",
               values_to = "Estimated_Flow") %>%
  mutate(Source_Site = case_when(
    Source_Site == "DAR_Q_08065800" ~ "Estimated flow using USGS-08065800",
    Source_Site == "DAR_Q_08109800" ~ "Estimated flow using USGS-08109800",
    Source_Site == "DAR_Q_08110100" ~ "Estimated flow using USGS-08110100"
  )) %>%
  ggplot() +
  geom_point(aes(Estimated_Flow, adjusted_flow, color = "DAR Estimates against Naturalized Flow"), alpha = 0.3) +
  geom_smooth(aes(Estimated_Flow, adjusted_flow, color = "DAR Estimates against Naturalized Flow"), 
              method = "lm", se = FALSE, alpha = 0.3) +
  geom_abline(aes(linetype = "1:1 Line", intercept = 0, slope = 1)) +
  facet_wrap(~Source_Site, scales = "free", ncol = 1) +
  scale_x_continuous(trans = "pseudo_log") + scale_y_continuous(trans = "pseudo_log") +
  labs(x = "DAR Estimated Flow [cfs]", y = "Naturalized Flow [cfs]") +
  guides(color = guide_legend(""), linetype = guide_legend("")) +
  theme_ms() + 
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "vertical",
        legend.box.margin = margin(0,0,0,0),
        legend.direction = "vertical") -> p1

dar_results_16396 %>%
  dplyr::select(Site, date, DAR_Q_08065800, DAR_Q_08109800, DAR_Q_08110100) %>%
  pivot_longer(c(DAR_Q_08065800, DAR_Q_08109800, DAR_Q_08110100),
               names_to = "Source_Site",
               values_to = "Estimated_Flow") %>%
  mutate(Source_Site = case_when(
    Source_Site == "DAR_Q_08065800" ~ "Estimated flow using USGS-08065800",
    Source_Site == "DAR_Q_08109800" ~ "Estimated flow using USGS-08109800",
    Source_Site == "DAR_Q_08110100" ~ "Estimated flow using USGS-08110100"
  )) %>%
  ggplot() +
  geom_line(aes(date, Estimated_Flow, color = "DAR Estimated Flow")) +
  geom_line(data = dar_results_16396, aes(date, adjusted_flow, color = "Naturalized Flow SWQM 16396"), alpha = 0.4) +
  facet_wrap(~Source_Site, scales = "free_y", ncol = 1) +
  labs(x = "Date", y = "Mean Daily Flow [cfs]") +
  theme_ms() + 
  theme(legend.title = element_blank(),
        legend.direction = "vertical") -> p2

p1 + p2
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<div class="figure">
<img src="document_files/figure-html/dar16396results-1.png" alt="Estimated streamflows at SWQM 16396 using DAR at three USGS gages plotted against measured streamflows and over time." width="100%" />
<p class="caption">(\#fig:dar16396results)Estimated streamflows at SWQM 16396 using DAR at three USGS gages plotted against measured streamflows and over time.</p>
</div>



```r
dar_results_16397 %>%
  dplyr::select(Site, date, adjusted_flow, DAR_Q_08065800, DAR_Q_08109800, DAR_Q_08110100) %>%
  pivot_longer(c(DAR_Q_08065800, DAR_Q_08109800, DAR_Q_08110100),
               names_to = "Source_Site",
               values_to = "Estimated_Flow") %>%
  mutate(Source_Site = case_when(
    Source_Site == "DAR_Q_08065800" ~ "Estimated flow using USGS-08065800",
    Source_Site == "DAR_Q_08109800" ~ "Estimated flow using USGS-08109800",
    Source_Site == "DAR_Q_08110100" ~ "Estimated flow using USGS-08110100"
  )) %>%
  ggplot() +
  geom_point(aes(Estimated_Flow, adjusted_flow, color = "DAR Estimates against Naturalized Flow"), alpha = 0.3) +
  geom_smooth(aes(Estimated_Flow, adjusted_flow, color = "DAR Estimates against Naturalized Flow"), 
              method = "lm", se = FALSE, alpha = 0.3) +
  geom_abline(aes(linetype = "1:1 Line", intercept = 0, slope = 1)) +
  facet_wrap(~Source_Site, scales = "free", ncol = 1) +
  scale_x_continuous(trans = "pseudo_log") + scale_y_continuous(trans = "pseudo_log") +
  labs(x = "DAR Estimated Flow [cfs]", y = "Naturalized Flow [cfs]") +
  guides(color = guide_legend(""), linetype = guide_legend("")) +
  theme_ms() + 
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "vertical",
        legend.box.margin = margin(0,0,0,0),
        legend.direction = "vertical") -> p1

dar_results_16397 %>%
  dplyr::select(Site, date, DAR_Q_08065800, DAR_Q_08109800, DAR_Q_08110100) %>%
  pivot_longer(c(DAR_Q_08065800, DAR_Q_08109800, DAR_Q_08110100),
               names_to = "Source_Site",
               values_to = "Estimated_Flow") %>%
  mutate(Source_Site = case_when(
    Source_Site == "DAR_Q_08065800" ~ "Estimated flow using USGS-08065800",
    Source_Site == "DAR_Q_08109800" ~ "Estimated flow using USGS-08109800",
    Source_Site == "DAR_Q_08110100" ~ "Estimated flow using USGS-08110100"
  )) %>%
  ggplot() +
  geom_line(aes(date, Estimated_Flow, color = "DAR Estimated Flow")) +
  geom_line(data = dar_results_16397, aes(date, adjusted_flow, color = "Naturalized Flow SWQM 16397"), alpha = 0.4) +
  facet_wrap(~Source_Site, scales = "free_y", ncol = 1) +
  labs(x = "Date", y = "Mean Daily Flow [cfs]") +
  theme_ms() + 
  theme(legend.title = element_blank(),
        legend.direction = "vertical") -> p2

p1 + p2
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<div class="figure">
<img src="document_files/figure-html/dar16397results-1.png" alt="Estimated streamflows at SWQM 16397 using DAR at three USGS gages plotted against measured streamflows and over time." width="100%" />
<p class="caption">(\#fig:dar16397results)Estimated streamflows at SWQM 16397 using DAR at three USGS gages plotted against measured streamflows and over time.</p>
</div>



```r
dar_results_16882 %>%
  dplyr::select(Site, date, adjusted_flow, DAR_Q_08065800, DAR_Q_08109800, DAR_Q_08110100) %>%
  pivot_longer(c(DAR_Q_08065800, DAR_Q_08109800, DAR_Q_08110100),
               names_to = "Source_Site",
               values_to = "Estimated_Flow") %>%
  mutate(Source_Site = case_when(
    Source_Site == "DAR_Q_08065800" ~ "Estimated flow using USGS-08065800",
    Source_Site == "DAR_Q_08109800" ~ "Estimated flow using USGS-08109800",
    Source_Site == "DAR_Q_08110100" ~ "Estimated flow using USGS-08110100"
  )) %>%
  ggplot() +
  geom_point(aes(Estimated_Flow, adjusted_flow, color = "DAR Estimates against Naturalized Flow"), alpha = 0.3) +
  geom_smooth(aes(Estimated_Flow, adjusted_flow, color = "DAR Estimates against Naturalized Flow"), 
              method = "lm", se = FALSE, alpha = 0.3) +
  geom_abline(aes(linetype = "1:1 Line", intercept = 0, slope = 1)) +
  facet_wrap(~Source_Site, scales = "free", ncol = 1) +
  scale_x_continuous(trans = "pseudo_log") + scale_y_continuous(trans = "pseudo_log") +
  labs(x = "DAR Estimated Flow [cfs]", y = "Naturalized Flow [cfs]") +
  guides(color = guide_legend(""), linetype = guide_legend("")) +
  theme_ms() + 
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "vertical",
        legend.box.margin = margin(0,0,0,0),
        legend.direction = "vertical") -> p1

dar_results_16882 %>%
  dplyr::select(Site, date, DAR_Q_08065800, DAR_Q_08109800, DAR_Q_08110100) %>%
  pivot_longer(c(DAR_Q_08065800, DAR_Q_08109800, DAR_Q_08110100),
               names_to = "Source_Site",
               values_to = "Estimated_Flow") %>%
  mutate(Source_Site = case_when(
    Source_Site == "DAR_Q_08065800" ~ "Estimated flow using USGS-08065800",
    Source_Site == "DAR_Q_08109800" ~ "Estimated flow using USGS-08109800",
    Source_Site == "DAR_Q_08110100" ~ "Estimated flow using USGS-08110100"
  )) %>%
  ggplot() +
  geom_line(aes(date, Estimated_Flow, color = "DAR Estimated Flow")) +
  geom_line(data = dar_results_16882, aes(date, adjusted_flow, color = "Naturalized Flow SWQM 16882"), alpha = 0.4) +
  facet_wrap(~Source_Site, scales = "free_y", ncol = 1) +
  labs(x = "Date", y = "Mean Daily Flow [cfs]") +
  theme_ms() + 
  theme(legend.title = element_blank(),
        legend.direction = "vertical") -> p2

p1 + p2
```

```
## `geom_smooth()` using formula 'y ~ x'
```

<div class="figure">
<img src="document_files/figure-html/dar16882results-1.png" alt="Estimated streamflows at SWQM 16882 using DAR at three USGS gages plotted against measured streamflows and over time." width="100%" />
<p class="caption">(\#fig:dar16882results)Estimated streamflows at SWQM 16882 using DAR at three USGS gages plotted against measured streamflows and over time.</p>
</div>

**Linear Regression**


```r
df_16396 %>%
  dplyr::select(Site, date, adjusted_flow) %>%
  left_join(usgs_08065800, by = c("date" = "datetime")) %>%
  left_join(usgs_08109800,  by = c("date" = "datetime")) %>%
  left_join(usgs_08110100,  by = c("date" = "datetime")) %>%
  mutate(lag_Flow_08065800 = lag(Flow_08065800),
         lag_Flow_08109800 = lag(Flow_08109800),
         lag_Flow_08110100 = lag(Flow_08110100),
         log_Q = log1p(adjusted_flow)) -> df_16396_lm

m1.lm <- lm(log_Q ~ Flow_08065800 + 
              Flow_08109800 + 
              Flow_08110100 + 
              lag_Flow_08065800 + 
              lag_Flow_08109800 + 
              lag_Flow_08110100,
            data = df_16396_lm)


flextable::as_flextable(m1.lm) %>%
  set_caption("Summary of linear regression coefficients at SWQM 16396")
```

```{=html}
<template id="9e8b8711-4ade-421c-a787-6f733f61b163"><style>
.tabwid table{
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 1.275em;
  margin-bottom: 1.275em;
  border-spacing: 0;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-b6b7329a{border-collapse:collapse;}.cl-b6a00480{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-b6a00481{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:italic;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-b6a0799c{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-b6a0799d{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-b6a13cce{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a13ccf{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a13cd0{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a13cd1{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a13cd2{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a13cd3{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a13cd4{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a13cd5{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a13cd6{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a13cd7{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a13cd8{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a163de{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a163df{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a163e0{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a163e1{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a163e2{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a163e3{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a163e4{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a163e5{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a163e6{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a163e7{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a163e8{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a18ada{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a18adb{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a18adc{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a18add{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a18ade{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a18adf{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a18ae0{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a18ae1{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a18ae2{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a18ae3{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a18ae4{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a1b1ea{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6a1b1eb{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-b6b7329a'>
```
<caption class="Table Caption">(\#tab:lrresults16396)Summary of linear regression coefficients at SWQM 16396</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-b6a1b1ea"><p class="cl-b6a0799c"><span class="cl-b6a00480"></span></p></td><td class="cl-b6a18ae3"><p class="cl-b6a0799d"><span class="cl-b6a00480">Estimate</span></p></td><td class="cl-b6a18ae4"><p class="cl-b6a0799d"><span class="cl-b6a00480">Standard Error</span></p></td><td class="cl-b6a1b1eb"><p class="cl-b6a0799d"><span class="cl-b6a00480">t value</span></p></td><td class="cl-b6a1b1eb"><p class="cl-b6a0799d"><span class="cl-b6a00480">Pr(&gt;|t|)</span></p></td><td class="cl-b6a18ae2"><p class="cl-b6a0799c"><span class="cl-b6a00480"></span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-b6a13cd2"><p class="cl-b6a0799c"><span class="cl-b6a00480">(Intercept)</span></p></td><td class="cl-b6a13cd1"><p class="cl-b6a0799d"><span class="cl-b6a00480">1.826</span></p></td><td class="cl-b6a13cd0"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.044</span></p></td><td class="cl-b6a13ccf"><p class="cl-b6a0799d"><span class="cl-b6a00480">41.191</span></p></td><td class="cl-b6a13ccf"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.0000</span></p></td><td class="cl-b6a13cce"><p class="cl-b6a0799c"><span class="cl-b6a00480">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6a13cd7"><p class="cl-b6a0799c"><span class="cl-b6a00480">Flow_08065800</span></p></td><td class="cl-b6a13cd3"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.000</span></p></td><td class="cl-b6a13cd4"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.000</span></p></td><td class="cl-b6a13cd5"><p class="cl-b6a0799d"><span class="cl-b6a00480">1.866</span></p></td><td class="cl-b6a13cd5"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.0631</span></p></td><td class="cl-b6a13cd6"><p class="cl-b6a0799c"><span class="cl-b6a00480">.</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6a13cd7"><p class="cl-b6a0799c"><span class="cl-b6a00480">Flow_08109800</span></p></td><td class="cl-b6a13cd3"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.017</span></p></td><td class="cl-b6a13cd4"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.002</span></p></td><td class="cl-b6a13cd5"><p class="cl-b6a0799d"><span class="cl-b6a00480">7.537</span></p></td><td class="cl-b6a13cd5"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.0000</span></p></td><td class="cl-b6a13cd6"><p class="cl-b6a0799c"><span class="cl-b6a00480">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6a13cd7"><p class="cl-b6a0799c"><span class="cl-b6a00480">Flow_08110100</span></p></td><td class="cl-b6a13cd3"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.005</span></p></td><td class="cl-b6a13cd4"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.001</span></p></td><td class="cl-b6a13cd5"><p class="cl-b6a0799d"><span class="cl-b6a00480">6.213</span></p></td><td class="cl-b6a13cd5"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.0000</span></p></td><td class="cl-b6a13cd6"><p class="cl-b6a0799c"><span class="cl-b6a00480">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6a163e1"><p class="cl-b6a0799c"><span class="cl-b6a00480">lag_Flow_08065800</span></p></td><td class="cl-b6a163e0"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.000</span></p></td><td class="cl-b6a163df"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.000</span></p></td><td class="cl-b6a163de"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.083</span></p></td><td class="cl-b6a163de"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.9342</span></p></td><td class="cl-b6a13cd8"><p class="cl-b6a0799c"><span class="cl-b6a00480"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6a163e1"><p class="cl-b6a0799c"><span class="cl-b6a00480">lag_Flow_08109800</span></p></td><td class="cl-b6a163e0"><p class="cl-b6a0799d"><span class="cl-b6a00480">-0.005</span></p></td><td class="cl-b6a163df"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.002</span></p></td><td class="cl-b6a163de"><p class="cl-b6a0799d"><span class="cl-b6a00480">-2.231</span></p></td><td class="cl-b6a163de"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.0264</span></p></td><td class="cl-b6a13cd8"><p class="cl-b6a0799c"><span class="cl-b6a00480">*</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6a163e4"><p class="cl-b6a0799c"><span class="cl-b6a00480">lag_Flow_08110100</span></p></td><td class="cl-b6a163e6"><p class="cl-b6a0799d"><span class="cl-b6a00480">-0.002</span></p></td><td class="cl-b6a163e5"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.001</span></p></td><td class="cl-b6a163e3"><p class="cl-b6a0799d"><span class="cl-b6a00480">-2.098</span></p></td><td class="cl-b6a163e3"><p class="cl-b6a0799d"><span class="cl-b6a00480">0.0367</span></p></td><td class="cl-b6a163e2"><p class="cl-b6a0799c"><span class="cl-b6a00480">*</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b6a18adc"><p class="cl-b6a0799d"><span class="cl-b6a00481">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b6a18adf"><p class="cl-b6a0799c"><span class="cl-b6a00480"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b6a18adf"><p class="cl-b6a0799c"><span class="cl-b6a00480">Residual standard error: 0.7083 on 304 degrees of freedom</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b6a18adf"><p class="cl-b6a0799c"><span class="cl-b6a00480">Multiple R-squared: 0.4088, Adjusted R-squared: 0.3972</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b6a18adf"><p class="cl-b6a0799c"><span class="cl-b6a00480">F-statistic: 35.04 on 304 and 6 DF, p-value: 0.0000</span></p></td></tr></tfoot></table></div></template>
<div class="flextable-shadow-host" id="18f94ea4-87e6-4a4f-8da6-f7debd780d5f"></div>
<script>
var dest = document.getElementById("18f94ea4-87e6-4a4f-8da6-f7debd780d5f");
var template = document.getElementById("9e8b8711-4ade-421c-a787-6f733f61b163");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;text-align:center;";
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

```


```r
df_16397 %>%
  dplyr::select(Site, date, adjusted_flow) %>%
  left_join(usgs_08065800, by = c("date" = "datetime")) %>%
  left_join(usgs_08109800,  by = c("date" = "datetime")) %>%
  left_join(usgs_08110100,  by = c("date" = "datetime")) %>%
  mutate(lag_Flow_08065800 = lag(Flow_08065800),
         lag_Flow_08109800 = lag(Flow_08109800),
         lag_Flow_08110100 = lag(Flow_08110100),
         log_Q = log1p(adjusted_flow)) -> df_16397_lm

m2.lm <- lm(log_Q ~ Flow_08065800 + 
              Flow_08109800 + 
              Flow_08110100 + 
              lag_Flow_08065800 + 
              lag_Flow_08109800 + 
              lag_Flow_08110100,
            data = df_16397_lm)


flextable::as_flextable(m2.lm) %>%
  set_caption("Summary of linear regression coefficients at SWQM 16397")
```

```{=html}
<template id="1af8c77f-f3e4-44b6-8825-d829916d87b7"><style>
.tabwid table{
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 1.275em;
  margin-bottom: 1.275em;
  border-spacing: 0;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-b75a388c{border-collapse:collapse;}.cl-b747c4a4{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-b747c4a5{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:italic;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-b74812ba{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-b74812bb{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-b74887e0{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b74887e1{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b74887e2{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b74887e3{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b74887e4{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b74887e5{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b74887e6{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b74887e7{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b74887e8{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b74887e9{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b74887ea{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748aee6{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748aee7{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748aee8{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748aee9{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748aeea{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748aeeb{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748aeec{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748aeed{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748aeee{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748aeef{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748aef0{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748d5e2{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748d5e3{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748d5e4{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748d5e5{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748d5e6{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748d5e7{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748d5e8{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748d5e9{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748d5ea{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748d5eb{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748d5ec{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748fcf2{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b748fcf3{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-b75a388c'>
```
<caption class="Table Caption">(\#tab:lrresults16397)Summary of linear regression coefficients at SWQM 16397</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-b748fcf2"><p class="cl-b74812ba"><span class="cl-b747c4a4"></span></p></td><td class="cl-b748d5eb"><p class="cl-b74812bb"><span class="cl-b747c4a4">Estimate</span></p></td><td class="cl-b748d5ec"><p class="cl-b74812bb"><span class="cl-b747c4a4">Standard Error</span></p></td><td class="cl-b748fcf3"><p class="cl-b74812bb"><span class="cl-b747c4a4">t value</span></p></td><td class="cl-b748fcf3"><p class="cl-b74812bb"><span class="cl-b747c4a4">Pr(&gt;|t|)</span></p></td><td class="cl-b748d5ea"><p class="cl-b74812ba"><span class="cl-b747c4a4"></span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-b74887e4"><p class="cl-b74812ba"><span class="cl-b747c4a4">(Intercept)</span></p></td><td class="cl-b74887e3"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.846</span></p></td><td class="cl-b74887e2"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.015</span></p></td><td class="cl-b74887e1"><p class="cl-b74812bb"><span class="cl-b747c4a4">57.074</span></p></td><td class="cl-b74887e1"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.0000</span></p></td><td class="cl-b74887e0"><p class="cl-b74812ba"><span class="cl-b747c4a4">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b74887e9"><p class="cl-b74812ba"><span class="cl-b747c4a4">Flow_08065800</span></p></td><td class="cl-b74887e5"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.000</span></p></td><td class="cl-b74887e6"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.000</span></p></td><td class="cl-b74887e7"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.695</span></p></td><td class="cl-b74887e7"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.4877</span></p></td><td class="cl-b74887e8"><p class="cl-b74812ba"><span class="cl-b747c4a4"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b74887e9"><p class="cl-b74812ba"><span class="cl-b747c4a4">Flow_08109800</span></p></td><td class="cl-b74887e5"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.006</span></p></td><td class="cl-b74887e6"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.001</span></p></td><td class="cl-b74887e7"><p class="cl-b74812bb"><span class="cl-b747c4a4">7.509</span></p></td><td class="cl-b74887e7"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.0000</span></p></td><td class="cl-b74887e8"><p class="cl-b74812ba"><span class="cl-b747c4a4">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b74887e9"><p class="cl-b74812ba"><span class="cl-b747c4a4">Flow_08110100</span></p></td><td class="cl-b74887e5"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.002</span></p></td><td class="cl-b74887e6"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.000</span></p></td><td class="cl-b74887e7"><p class="cl-b74812bb"><span class="cl-b747c4a4">8.296</span></p></td><td class="cl-b74887e7"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.0000</span></p></td><td class="cl-b74887e8"><p class="cl-b74812ba"><span class="cl-b747c4a4">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b748aee9"><p class="cl-b74812ba"><span class="cl-b747c4a4">lag_Flow_08065800</span></p></td><td class="cl-b748aee8"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.000</span></p></td><td class="cl-b748aee7"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.000</span></p></td><td class="cl-b748aee6"><p class="cl-b74812bb"><span class="cl-b747c4a4">1.739</span></p></td><td class="cl-b748aee6"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.0830</span></p></td><td class="cl-b74887ea"><p class="cl-b74812ba"><span class="cl-b747c4a4">.</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b748aee9"><p class="cl-b74812ba"><span class="cl-b747c4a4">lag_Flow_08109800</span></p></td><td class="cl-b748aee8"><p class="cl-b74812bb"><span class="cl-b747c4a4">-0.003</span></p></td><td class="cl-b748aee7"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.001</span></p></td><td class="cl-b748aee6"><p class="cl-b74812bb"><span class="cl-b747c4a4">-3.970</span></p></td><td class="cl-b748aee6"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.0001</span></p></td><td class="cl-b74887ea"><p class="cl-b74812ba"><span class="cl-b747c4a4">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b748aeec"><p class="cl-b74812ba"><span class="cl-b747c4a4">lag_Flow_08110100</span></p></td><td class="cl-b748aeee"><p class="cl-b74812bb"><span class="cl-b747c4a4">-0.001</span></p></td><td class="cl-b748aeed"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.000</span></p></td><td class="cl-b748aeeb"><p class="cl-b74812bb"><span class="cl-b747c4a4">-2.920</span></p></td><td class="cl-b748aeeb"><p class="cl-b74812bb"><span class="cl-b747c4a4">0.0038</span></p></td><td class="cl-b748aeea"><p class="cl-b74812ba"><span class="cl-b747c4a4">**</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b748d5e4"><p class="cl-b74812bb"><span class="cl-b747c4a5">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b748d5e7"><p class="cl-b74812ba"><span class="cl-b747c4a4"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b748d5e7"><p class="cl-b74812ba"><span class="cl-b747c4a4">Residual standard error: 0.2368 on 304 degrees of freedom</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b748d5e7"><p class="cl-b74812ba"><span class="cl-b747c4a4">Multiple R-squared: 0.4337, Adjusted R-squared: 0.4226</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b748d5e7"><p class="cl-b74812ba"><span class="cl-b747c4a4">F-statistic: 38.81 on 304 and 6 DF, p-value: 0.0000</span></p></td></tr></tfoot></table></div></template>
<div class="flextable-shadow-host" id="b298e56d-f6b1-46b4-9210-54ee8bc3312f"></div>
<script>
var dest = document.getElementById("b298e56d-f6b1-46b4-9210-54ee8bc3312f");
var template = document.getElementById("1af8c77f-f3e4-44b6-8825-d829916d87b7");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;text-align:center;";
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

```


```r
df_16882 %>%
  dplyr::select(Site, date, adjusted_flow) %>%
  left_join(usgs_08065800, by = c("date" = "datetime")) %>%
  left_join(usgs_08109800,  by = c("date" = "datetime")) %>%
  left_join(usgs_08110100,  by = c("date" = "datetime")) %>%
  mutate(lag_Flow_08065800 = lag(Flow_08065800),
         lag_Flow_08109800 = lag(Flow_08109800),
         lag_Flow_08110100 = lag(Flow_08110100),
         log_Q = log1p(adjusted_flow)) -> df_16882_lm

m3.lm <- lm(log_Q ~ Flow_08065800 + 
              Flow_08109800 + 
              Flow_08110100 + 
              lag_Flow_08065800 + 
              lag_Flow_08109800 + 
              lag_Flow_08110100,
            data = df_16882_lm)


flextable::as_flextable(m3.lm) %>%
  set_caption("Summary of linear regression coefficients at SWQM 16882")
```

```{=html}
<template id="adf83bbb-41d0-45bf-b183-adb07574b057"><style>
.tabwid table{
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 1.275em;
  margin-bottom: 1.275em;
  border-spacing: 0;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-b7f18070{border-collapse:collapse;}.cl-b7dff6b6{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-b7dff6b7{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:italic;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-b7e01dbc{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-b7e01dbd{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-b7e06bdc{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e06bdd{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e06bde{width:52.9pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e06bdf{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e06be0{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e06be1{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e06be2{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e06be3{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e06be4{width:52.9pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e06be5{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e06be6{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e092d8{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e092d9{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e092da{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e092db{width:52.9pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e092dc{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e092dd{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e092de{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e092df{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e092e0{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e092e1{width:52.9pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e092e2{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0b9e8{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0b9e9{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0b9ea{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0b9eb{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0b9ec{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0b9ed{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0b9ee{width:52.9pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0b9ef{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0b9f0{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0b9f1{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0b9f2{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0e0ee{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0e0ef{width:52.9pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0e0f0{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0e0f1{width:33.4pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0e0f2{width:63.3pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0e0f3{width:92.7pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0e0f4{width:119.6pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0e0f5{width:52.9pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b7e0e0f6{width:54.2pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(102, 102, 102, 1.00);border-top: 2pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-b7f18070'>
```
<caption class="Table Caption">(\#tab:lrresults16882)Summary of linear regression coefficients at SWQM 16882</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-b7e0e0f4"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6"></span></p></td><td class="cl-b7e0e0f2"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">Estimate</span></p></td><td class="cl-b7e0e0f3"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">Standard Error</span></p></td><td class="cl-b7e0e0f5"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">t value</span></p></td><td class="cl-b7e0e0f6"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">Pr(&gt;|t|)</span></p></td><td class="cl-b7e0e0f1"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6"></span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-b7e06be1"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">(Intercept)</span></p></td><td class="cl-b7e06be0"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.312</span></p></td><td class="cl-b7e06bdf"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.049</span></p></td><td class="cl-b7e06bde"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">6.354</span></p></td><td class="cl-b7e06bdd"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.0000</span></p></td><td class="cl-b7e06bdc"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b7e092d8"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">Flow_08065800</span></p></td><td class="cl-b7e06be2"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">-0.000</span></p></td><td class="cl-b7e06be3"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.000</span></p></td><td class="cl-b7e06be4"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">-0.087</span></p></td><td class="cl-b7e06be5"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.9309</span></p></td><td class="cl-b7e06be6"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b7e092d8"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">Flow_08109800</span></p></td><td class="cl-b7e06be2"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.018</span></p></td><td class="cl-b7e06be3"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.003</span></p></td><td class="cl-b7e06be4"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">7.185</span></p></td><td class="cl-b7e06be5"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.0000</span></p></td><td class="cl-b7e06be6"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b7e092d8"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">Flow_08110100</span></p></td><td class="cl-b7e06be2"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.004</span></p></td><td class="cl-b7e06be3"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.001</span></p></td><td class="cl-b7e06be4"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">3.933</span></p></td><td class="cl-b7e06be5"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.0001</span></p></td><td class="cl-b7e06be6"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b7e092de"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">lag_Flow_08065800</span></p></td><td class="cl-b7e092dd"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.000</span></p></td><td class="cl-b7e092dc"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.000</span></p></td><td class="cl-b7e092db"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.640</span></p></td><td class="cl-b7e092da"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.5223</span></p></td><td class="cl-b7e092d9"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b7e092de"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">lag_Flow_08109800</span></p></td><td class="cl-b7e092dd"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">-0.009</span></p></td><td class="cl-b7e092dc"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.002</span></p></td><td class="cl-b7e092db"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">-3.524</span></p></td><td class="cl-b7e092da"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.0005</span></p></td><td class="cl-b7e092d9"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b7e092e2"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">lag_Flow_08110100</span></p></td><td class="cl-b7e0b9e9"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">-0.003</span></p></td><td class="cl-b7e0b9e8"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.001</span></p></td><td class="cl-b7e092e1"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">-2.702</span></p></td><td class="cl-b7e092e0"><p class="cl-b7e01dbd"><span class="cl-b7dff6b6">0.0073</span></p></td><td class="cl-b7e092df"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">**</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b7e0b9ef"><p class="cl-b7e01dbd"><span class="cl-b7dff6b7">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b7e0b9f2"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b7e0b9f2"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">Residual standard error: 0.7923 on 310 degrees of freedom</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b7e0b9f2"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">Multiple R-squared: 0.2595, Adjusted R-squared: 0.2452</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b7e0b9f2"><p class="cl-b7e01dbc"><span class="cl-b7dff6b6">F-statistic: 18.11 on 310 and 6 DF, p-value: 0.0000</span></p></td></tr></tfoot></table></div></template>
<div class="flextable-shadow-host" id="b96cd84d-38e0-4995-8df2-3e33eb1a6612"></div>
<script>
var dest = document.getElementById("b96cd84d-38e0-4995-8df2-3e33eb1a6612");
var template = document.getElementById("adf83bbb-41d0-45bf-b183-adb07574b057");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;text-align:center;";
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

```


```r
df_16396_lm %>%
  mutate(fits = as.numeric(predict(m1.lm, newdata = .))) %>%
  mutate(fits = expm1(fits)) -> df_16396_lm


df_16397_lm %>%
  mutate(fits = as.numeric(predict(m2.lm, newdata = .))) %>%
  mutate(fits = expm1(fits)) -> df_16397_lm

df_16882_lm %>%
  mutate(fits = as.numeric(predict(m3.lm, newdata = .))) %>%
  mutate(fits = expm1(fits)) -> df_16882_lm

df_16396_lm %>%
  ggplot() +
  geom_point(aes(fits, adjusted_flow, color = "Regression Estimates against Naturalized Flow"),
             alpha =0.4) +
  geom_smooth(aes(fits, adjusted_flow, color = "Regression Estimates against Naturalized Flow"),
              method = "lm", se = FALSE) +
  geom_abline(aes(linetype = "1:1 Line", intercept = 0, slope = 1)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Regression Estimated Flow [cfs]", y = "Naturalized Flow [cfs]",
       title = "SWQM-16396") +
  guides(color = guide_legend(""), linetype = guide_legend("")) +
  coord_equal() +
  theme_ms() +
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "vertical",
        legend.box.margin = margin(0,0,0,0),
        legend.direction = "vertical",
        legend.position = "none")-> p1

df_16396_lm %>%
  ggplot() +
  geom_line(aes(date, fits, color = "Linear Regression Estimated Flow")) +
  geom_line(aes(date, adjusted_flow, color = "Naturalized Flow SWQM 16396"), alpha = 0.4) +
  labs(x = "Date", y = "Mean Daily Flow [cfs]") +
  theme_ms() +
  theme(legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "none") -> p2


df_16397_lm %>%
  ggplot() +
  geom_point(aes(fits, adjusted_flow, color = "Regression Estimates against Naturalized Flow"),
             alpha =0.4) +
  geom_smooth(aes(fits, adjusted_flow, color = "Regression Estimates against Naturalized Flow"),
              method = "lm", se = FALSE) +
  geom_abline(aes(linetype = "1:1 Line", intercept = 0, slope = 1)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Regression Estimated Flow [cfs]", y = "Naturalized Flow [cfs]",
       title = "SWQM-16397") +
  guides(color = guide_legend(""), linetype = guide_legend("")) +
  coord_equal() +
  theme_ms() +
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "vertical",
        legend.box.margin = margin(0,0,0,0),
        legend.direction = "vertical",
        legend.position = "none")-> p3

df_16397_lm %>%
  ggplot() +
  geom_line(aes(date, fits, color = "Linear Regression Estimated Flow")) +
  geom_line(aes(date, adjusted_flow, color = "Naturalized Flow SWQM 16397"), alpha = 0.4) +
  labs(x = "Date", y = "Mean Daily Flow [cfs]") +
  theme_ms() +
  theme(legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "none") -> p4

df_16882_lm %>%
  ggplot() +
  geom_point(aes(fits, adjusted_flow, color = "Regression Estimates against Naturalized Flow"),
             alpha =0.4) +
  geom_smooth(aes(fits, adjusted_flow, color = "Regression Estimates against Naturalized Flow"),
              method = "lm", se = FALSE) +
  geom_abline(aes(linetype = "1:1 Line", intercept = 0, slope = 1)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Regression Estimated Flow [cfs]", y = "Naturalized Flow [cfs]",
       title = "SWQM-16882") +
  guides(color = guide_legend(""), linetype = guide_legend("")) +
  coord_equal() +
  theme_ms() +
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "vertical",
        legend.box.margin = margin(0,0,0,0),
        legend.direction = "vertical")-> p5

df_16882_lm %>%
  ggplot() +
  geom_line(aes(date, fits, color = "Linear Regression Estimated Flow")) +
  geom_line(aes(date, adjusted_flow, color = "Naturalized Flow"), alpha = 0.4) +
  labs(x = "Date", y = "Mean Daily Flow [cfs]") +
  theme_ms() +
  theme(legend.title = element_blank(),
        legend.direction = "vertical") -> p6

(p1 + p2)/(p3 + p4)/(p5 + p6)
```

<div class="figure">
<img src="document_files/figure-html/unnamed-chunk-4-1.png" alt="Linear regression estimated mean daily streamflows plotted against measured mean daily streamflows and over time." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-4)Linear regression estimated mean daily streamflows plotted against measured mean daily streamflows and over time.</p>
</div>


**GAM**


```r
m1.gam <- gam(adjusted_flow ~
            s(ewood_precip, bs = "ts", k = 5) +
            s(ewood_tmax, bs = "ts", k = 5) +
            s(lagPrecip, bs = "ts", k = 5) +
            s(doy, bs = "cc", k = 8) +
            s(wetness, bs = "ts", k = 5) +
            s(et, bs = "ts", k = 5) +
              ti(ewood_precip, lagPrecip, k = 5) +
              ti(ewood_precip, ewood_tmax, k = 5) +
              ti(ewood_precip, wetness, k = 5) +
              ti(ewood_precip, et, k = 5) +
              ti(ewood_tmax, wetness, k = 5) +
              ti(ewood_tmax, et, k = 5) +
              ti(wetness, et, k = 5),
          data = df_16396,
          select = TRUE,
          family = scat(link = "log"),
          method = "REML",
          control = gam.control(nthreads = 2))
flextable::as_flextable(m1.gam) %>%
  flextable::set_caption("Summary of GAM coefficients at SWQM 16396")
```

```{=html}
<template id="49bdf7a6-21f0-4951-854c-74fd46dbfe5d"><style>
.tabwid table{
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 1.275em;
  margin-bottom: 1.275em;
  border-spacing: 0;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-4a3db210{border-collapse:collapse;}.cl-4a299a64{font-family:'Arial';font-size:11pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-4a299a65{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-4a299a66{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4a299a67{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4a299a68{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4a299a69{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a299a6a{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a299a6b{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a299a6c{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a299a6d{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a299a6e{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2bfc8c{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2bfc8d{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2bfc8e{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2bfc8f{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2bfc90{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2bfc91{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2bfc92{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2bfc93{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2bfc94{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2bfc95{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2bfc96{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2e6396{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2e6397{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2e6398{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2e6399{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2e639a{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2e639b{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2e639c{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2e639d{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2e639e{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2e639f{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2e63a0{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2fea68{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2fea69{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2fea6a{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2fea6b{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2fea6c{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2fea6d{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2fea6e{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2fea6f{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2fea70{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2fea71{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a2fea72{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a324da8{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a324da9{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4a324daa{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-4a3db210'>
```
<caption class="Table Caption">(\#tab:gamresults16396)Summary of GAM coefficients at SWQM 16396</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-4a2fea70"><p class="cl-4a299a66"><span class="cl-4a299a64">Component</span></p></td><td class="cl-4a324da8"><p class="cl-4a299a66"><span class="cl-4a299a64">Term</span></p></td><td class="cl-4a2fea71"><p class="cl-4a299a67"><span class="cl-4a299a64">Estimate</span></p></td><td class="cl-4a324daa"><p class="cl-4a299a67"><span class="cl-4a299a64">Std Error</span></p></td><td class="cl-4a324da9"><p class="cl-4a299a67"><span class="cl-4a299a64">t-value</span></p></td><td class="cl-4a2fea72"><p class="cl-4a299a67"><span class="cl-4a299a64">p-value</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-4a299a6d"><p class="cl-4a299a68"><span class="cl-4a299a65">A. parametric coefficients</span></p></td><td class="cl-4a299a6c"><p class="cl-4a299a66"><span class="cl-4a299a65">(Intercept)</span></p></td><td class="cl-4a299a6b"><p class="cl-4a299a67"><span class="cl-4a299a65">1.879</span></p></td><td class="cl-4a299a6a"><p class="cl-4a299a67"><span class="cl-4a299a65">0.052</span></p></td><td class="cl-4a299a69"><p class="cl-4a299a67"><span class="cl-4a299a65">36.220</span></p></td><td class="cl-4a299a6e"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2e639a"><p class="cl-4a299a68"><span class="cl-4a299a64">Component</span></p></td><td class="cl-4a2e6398"><p class="cl-4a299a66"><span class="cl-4a299a64">Term</span></p></td><td class="cl-4a2e6399"><p class="cl-4a299a67"><span class="cl-4a299a64">edf</span></p></td><td class="cl-4a2e6397"><p class="cl-4a299a67"><span class="cl-4a299a64">Ref. df</span></p></td><td class="cl-4a2e639b"><p class="cl-4a299a67"><span class="cl-4a299a64">F-value</span></p></td><td class="cl-4a2e639c"><p class="cl-4a299a67"><span class="cl-4a299a64">p-value</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  rowspan="13"class="cl-4a2fea69"><p class="cl-4a299a68"><span class="cl-4a299a65">B. smooth terms</span></p></td><td class="cl-4a2fea68"><p class="cl-4a299a66"><span class="cl-4a299a65">s(ewood_precip)</span></p></td><td class="cl-4a2e639e"><p class="cl-4a299a67"><span class="cl-4a299a65">0.919</span></p></td><td class="cl-4a2e639d"><p class="cl-4a299a67"><span class="cl-4a299a65">4.000</span></p></td><td class="cl-4a2e63a0"><p class="cl-4a299a67"><span class="cl-4a299a65">9.532</span></p></td><td class="cl-4a2e639f"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2bfc8d"><p class="cl-4a299a66"><span class="cl-4a299a65">s(ewood_tmax)</span></p></td><td class="cl-4a2bfc91"><p class="cl-4a299a67"><span class="cl-4a299a65">0.984</span></p></td><td class="cl-4a2bfc8e"><p class="cl-4a299a67"><span class="cl-4a299a65">4.000</span></p></td><td class="cl-4a2bfc8f"><p class="cl-4a299a67"><span class="cl-4a299a65">9.149</span></p></td><td class="cl-4a2bfc8c"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2bfc8d"><p class="cl-4a299a66"><span class="cl-4a299a65">s(lagPrecip)</span></p></td><td class="cl-4a2bfc91"><p class="cl-4a299a67"><span class="cl-4a299a65">1.140</span></p></td><td class="cl-4a2bfc8e"><p class="cl-4a299a67"><span class="cl-4a299a65">4.000</span></p></td><td class="cl-4a2bfc8f"><p class="cl-4a299a67"><span class="cl-4a299a65">309.675</span></p></td><td class="cl-4a2bfc8c"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2bfc8d"><p class="cl-4a299a66"><span class="cl-4a299a65">s(doy)</span></p></td><td class="cl-4a2bfc91"><p class="cl-4a299a67"><span class="cl-4a299a65">5.138</span></p></td><td class="cl-4a2bfc8e"><p class="cl-4a299a67"><span class="cl-4a299a65">6.000</span></p></td><td class="cl-4a2bfc8f"><p class="cl-4a299a67"><span class="cl-4a299a65">333.401</span></p></td><td class="cl-4a2bfc8c"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2bfc8d"><p class="cl-4a299a66"><span class="cl-4a299a65">s(wetness)</span></p></td><td class="cl-4a2bfc91"><p class="cl-4a299a67"><span class="cl-4a299a65">1.023</span></p></td><td class="cl-4a2bfc8e"><p class="cl-4a299a67"><span class="cl-4a299a65">4.000</span></p></td><td class="cl-4a2bfc8f"><p class="cl-4a299a67"><span class="cl-4a299a65">48.803</span></p></td><td class="cl-4a2bfc8c"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2bfc8d"><p class="cl-4a299a66"><span class="cl-4a299a65">s(et)</span></p></td><td class="cl-4a2bfc91"><p class="cl-4a299a67"><span class="cl-4a299a65">0.651</span></p></td><td class="cl-4a2bfc8e"><p class="cl-4a299a67"><span class="cl-4a299a65">4.000</span></p></td><td class="cl-4a2bfc8f"><p class="cl-4a299a67"><span class="cl-4a299a65">2.135</span></p></td><td class="cl-4a2bfc8c"><p class="cl-4a299a67"><span class="cl-4a299a65">*</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2bfc8d"><p class="cl-4a299a66"><span class="cl-4a299a65">ti(ewood_precip,lagPrecip)</span></p></td><td class="cl-4a2bfc91"><p class="cl-4a299a67"><span class="cl-4a299a65">5.144</span></p></td><td class="cl-4a2bfc8e"><p class="cl-4a299a67"><span class="cl-4a299a65">16.000</span></p></td><td class="cl-4a2bfc8f"><p class="cl-4a299a67"><span class="cl-4a299a65">18.487</span></p></td><td class="cl-4a2bfc8c"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2bfc8d"><p class="cl-4a299a66"><span class="cl-4a299a65">ti(ewood_precip,ewood_tmax)</span></p></td><td class="cl-4a2bfc91"><p class="cl-4a299a67"><span class="cl-4a299a65">9.684</span></p></td><td class="cl-4a2bfc8e"><p class="cl-4a299a67"><span class="cl-4a299a65">16.000</span></p></td><td class="cl-4a2bfc8f"><p class="cl-4a299a67"><span class="cl-4a299a65">130.341</span></p></td><td class="cl-4a2bfc8c"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2bfc8d"><p class="cl-4a299a66"><span class="cl-4a299a65">ti(ewood_precip,wetness)</span></p></td><td class="cl-4a2bfc91"><p class="cl-4a299a67"><span class="cl-4a299a65">8.479</span></p></td><td class="cl-4a2bfc8e"><p class="cl-4a299a67"><span class="cl-4a299a65">16.000</span></p></td><td class="cl-4a2bfc8f"><p class="cl-4a299a67"><span class="cl-4a299a65">133.633</span></p></td><td class="cl-4a2bfc8c"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2bfc8d"><p class="cl-4a299a66"><span class="cl-4a299a65">ti(ewood_precip,et)</span></p></td><td class="cl-4a2bfc91"><p class="cl-4a299a67"><span class="cl-4a299a65">11.560</span></p></td><td class="cl-4a2bfc8e"><p class="cl-4a299a67"><span class="cl-4a299a65">16.000</span></p></td><td class="cl-4a2bfc8f"><p class="cl-4a299a67"><span class="cl-4a299a65">773.912</span></p></td><td class="cl-4a2bfc8c"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2bfc8d"><p class="cl-4a299a66"><span class="cl-4a299a65">ti(ewood_tmax,wetness)</span></p></td><td class="cl-4a2bfc91"><p class="cl-4a299a67"><span class="cl-4a299a65">4.628</span></p></td><td class="cl-4a2bfc8e"><p class="cl-4a299a67"><span class="cl-4a299a65">16.000</span></p></td><td class="cl-4a2bfc8f"><p class="cl-4a299a67"><span class="cl-4a299a65">31.182</span></p></td><td class="cl-4a2bfc8c"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2bfc8d"><p class="cl-4a299a66"><span class="cl-4a299a65">ti(ewood_tmax,et)</span></p></td><td class="cl-4a2bfc91"><p class="cl-4a299a67"><span class="cl-4a299a65">6.571</span></p></td><td class="cl-4a2bfc8e"><p class="cl-4a299a67"><span class="cl-4a299a65">16.000</span></p></td><td class="cl-4a2bfc8f"><p class="cl-4a299a67"><span class="cl-4a299a65">39.683</span></p></td><td class="cl-4a2bfc8c"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4a2bfc94"><p class="cl-4a299a66"><span class="cl-4a299a65">ti(wetness,et)</span></p></td><td class="cl-4a2e6396"><p class="cl-4a299a67"><span class="cl-4a299a65">11.358</span></p></td><td class="cl-4a2bfc93"><p class="cl-4a299a67"><span class="cl-4a299a65">16.000</span></p></td><td class="cl-4a2bfc92"><p class="cl-4a299a67"><span class="cl-4a299a65">457.621</span></p></td><td class="cl-4a2bfc95"><p class="cl-4a299a67"><span class="cl-4a299a65">***</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4a2fea6f"><p class="cl-4a299a66"><span class="cl-4a299a65">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4a2fea6f"><p class="cl-4a299a66"><span class="cl-4a299a65"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4a2fea6f"><p class="cl-4a299a66"><span class="cl-4a299a65">Adjusted R-squared: 0.494, Deviance explained 0.632</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4a2fea6f"><p class="cl-4a299a66"><span class="cl-4a299a65">-REML : 890.024, Scale est: 1.000, N: 312</span></p></td></tr></tfoot></table></div></template>
<div class="flextable-shadow-host" id="b77729bd-b43c-4a91-a7a6-c24959a4f540"></div>
<script>
var dest = document.getElementById("b77729bd-b43c-4a91-a7a6-c24959a4f540");
var template = document.getElementById("49bdf7a6-21f0-4951-854c-74fd46dbfe5d");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;text-align:center;";
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

```



```r
m2.gam <- gam(adjusted_flow ~
            s(ewood_precip, bs = "ts", k = 5) +
            s(ewood_tmax, bs = "ts", k = 5) +
            s(lagPrecip, bs = "ts", k = 5) +
            s(doy, bs = "cc", k = 8) +
            s(wetness, bs = "ts", k = 5) +
            s(et, bs = "ts", k = 5) +
              ti(ewood_precip, lagPrecip, k = 5) +
              ti(ewood_precip, ewood_tmax, k = 5) +
              ti(ewood_precip, wetness, k = 5) +
              ti(ewood_precip, et, k = 5) +
              ti(ewood_tmax, wetness, k = 5) +
              ti(ewood_tmax, et, k = 5) +
              ti(wetness, et, k = 5),
          data = df_16397,
          select = TRUE,
          family = scat(link = "log"),
          method = "REML",
          control = gam.control(nthreads = 2))
flextable::as_flextable(m2.gam) %>%
  flextable::set_caption("Summary of GAM coefficients at SWQM 16397")
```

```{=html}
<template id="1e1f21ef-add5-4b6d-a039-3379ea594b8c"><style>
.tabwid table{
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 1.275em;
  margin-bottom: 1.275em;
  border-spacing: 0;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-d253ba2c{border-collapse:collapse;}.cl-d23fa78a{font-family:'Arial';font-size:11pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-d23fa78b{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-d23fa78c{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-d23fa78d{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-d23fa78e{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-d23fa78f{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d23fa790{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d23fa791{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d23fa792{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d23fa793{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d23fa794{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d24209da{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d24209db{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d24209dc{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d24209dd{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d24209de{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d24209df{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d24209e0{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d24209e1{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d24209e2{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d24209e3{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d24209e4{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2446c0c{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2446c0d{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2446c0e{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2446c0f{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2446c10{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2446c11{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2446c12{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2446c13{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2446c14{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2446c15{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2446c16{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d246ce52{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d246ce53{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d246ce54{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d246ce55{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d246ce56{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d246ce57{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d246ce58{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d246ce59{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d246ce5a{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d246ce5b{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d246ce5c{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2492fe4{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2492fe5{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-d2492fe6{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-d253ba2c'>
```
<caption class="Table Caption">(\#tab:gamresults16397)Summary of GAM coefficients at SWQM 16397</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-d246ce5a"><p class="cl-d23fa78c"><span class="cl-d23fa78a">Component</span></p></td><td class="cl-d2492fe4"><p class="cl-d23fa78c"><span class="cl-d23fa78a">Term</span></p></td><td class="cl-d246ce5b"><p class="cl-d23fa78d"><span class="cl-d23fa78a">Estimate</span></p></td><td class="cl-d2492fe6"><p class="cl-d23fa78d"><span class="cl-d23fa78a">Std Error</span></p></td><td class="cl-d2492fe5"><p class="cl-d23fa78d"><span class="cl-d23fa78a">t-value</span></p></td><td class="cl-d246ce5c"><p class="cl-d23fa78d"><span class="cl-d23fa78a">p-value</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-d23fa793"><p class="cl-d23fa78e"><span class="cl-d23fa78b">A. parametric coefficients</span></p></td><td class="cl-d23fa792"><p class="cl-d23fa78c"><span class="cl-d23fa78b">(Intercept)</span></p></td><td class="cl-d23fa791"><p class="cl-d23fa78d"><span class="cl-d23fa78b">0.405</span></p></td><td class="cl-d23fa790"><p class="cl-d23fa78d"><span class="cl-d23fa78b">0.012</span></p></td><td class="cl-d23fa78f"><p class="cl-d23fa78d"><span class="cl-d23fa78b">32.636</span></p></td><td class="cl-d23fa794"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d2446c10"><p class="cl-d23fa78e"><span class="cl-d23fa78a">Component</span></p></td><td class="cl-d2446c0e"><p class="cl-d23fa78c"><span class="cl-d23fa78a">Term</span></p></td><td class="cl-d2446c0f"><p class="cl-d23fa78d"><span class="cl-d23fa78a">edf</span></p></td><td class="cl-d2446c0d"><p class="cl-d23fa78d"><span class="cl-d23fa78a">Ref. df</span></p></td><td class="cl-d2446c11"><p class="cl-d23fa78d"><span class="cl-d23fa78a">F-value</span></p></td><td class="cl-d2446c12"><p class="cl-d23fa78d"><span class="cl-d23fa78a">p-value</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  rowspan="13"class="cl-d246ce53"><p class="cl-d23fa78e"><span class="cl-d23fa78b">B. smooth terms</span></p></td><td class="cl-d246ce52"><p class="cl-d23fa78c"><span class="cl-d23fa78b">s(ewood_precip)</span></p></td><td class="cl-d2446c14"><p class="cl-d23fa78d"><span class="cl-d23fa78b">0.917</span></p></td><td class="cl-d2446c13"><p class="cl-d23fa78d"><span class="cl-d23fa78b">4.000</span></p></td><td class="cl-d2446c16"><p class="cl-d23fa78d"><span class="cl-d23fa78b">12.252</span></p></td><td class="cl-d2446c15"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d24209db"><p class="cl-d23fa78c"><span class="cl-d23fa78b">s(ewood_tmax)</span></p></td><td class="cl-d24209df"><p class="cl-d23fa78d"><span class="cl-d23fa78b">0.683</span></p></td><td class="cl-d24209dc"><p class="cl-d23fa78d"><span class="cl-d23fa78b">4.000</span></p></td><td class="cl-d24209dd"><p class="cl-d23fa78d"><span class="cl-d23fa78b">1.546</span></p></td><td class="cl-d24209da"><p class="cl-d23fa78d"><span class="cl-d23fa78b">*</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d24209db"><p class="cl-d23fa78c"><span class="cl-d23fa78b">s(lagPrecip)</span></p></td><td class="cl-d24209df"><p class="cl-d23fa78d"><span class="cl-d23fa78b">1.033</span></p></td><td class="cl-d24209dc"><p class="cl-d23fa78d"><span class="cl-d23fa78b">4.000</span></p></td><td class="cl-d24209dd"><p class="cl-d23fa78d"><span class="cl-d23fa78b">57.638</span></p></td><td class="cl-d24209da"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d24209db"><p class="cl-d23fa78c"><span class="cl-d23fa78b">s(doy)</span></p></td><td class="cl-d24209df"><p class="cl-d23fa78d"><span class="cl-d23fa78b">5.852</span></p></td><td class="cl-d24209dc"><p class="cl-d23fa78d"><span class="cl-d23fa78b">6.000</span></p></td><td class="cl-d24209dd"><p class="cl-d23fa78d"><span class="cl-d23fa78b">302.695</span></p></td><td class="cl-d24209da"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d24209db"><p class="cl-d23fa78c"><span class="cl-d23fa78b">s(wetness)</span></p></td><td class="cl-d24209df"><p class="cl-d23fa78d"><span class="cl-d23fa78b">1.029</span></p></td><td class="cl-d24209dc"><p class="cl-d23fa78d"><span class="cl-d23fa78b">4.000</span></p></td><td class="cl-d24209dd"><p class="cl-d23fa78d"><span class="cl-d23fa78b">33.736</span></p></td><td class="cl-d24209da"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d24209db"><p class="cl-d23fa78c"><span class="cl-d23fa78b">s(et)</span></p></td><td class="cl-d24209df"><p class="cl-d23fa78d"><span class="cl-d23fa78b">1.353</span></p></td><td class="cl-d24209dc"><p class="cl-d23fa78d"><span class="cl-d23fa78b">4.000</span></p></td><td class="cl-d24209dd"><p class="cl-d23fa78d"><span class="cl-d23fa78b">27.929</span></p></td><td class="cl-d24209da"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d24209db"><p class="cl-d23fa78c"><span class="cl-d23fa78b">ti(ewood_precip,lagPrecip)</span></p></td><td class="cl-d24209df"><p class="cl-d23fa78d"><span class="cl-d23fa78b">10.805</span></p></td><td class="cl-d24209dc"><p class="cl-d23fa78d"><span class="cl-d23fa78b">16.000</span></p></td><td class="cl-d24209dd"><p class="cl-d23fa78d"><span class="cl-d23fa78b">216.589</span></p></td><td class="cl-d24209da"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d24209db"><p class="cl-d23fa78c"><span class="cl-d23fa78b">ti(ewood_precip,ewood_tmax)</span></p></td><td class="cl-d24209df"><p class="cl-d23fa78d"><span class="cl-d23fa78b">10.582</span></p></td><td class="cl-d24209dc"><p class="cl-d23fa78d"><span class="cl-d23fa78b">16.000</span></p></td><td class="cl-d24209dd"><p class="cl-d23fa78d"><span class="cl-d23fa78b">168.673</span></p></td><td class="cl-d24209da"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d24209db"><p class="cl-d23fa78c"><span class="cl-d23fa78b">ti(ewood_precip,wetness)</span></p></td><td class="cl-d24209df"><p class="cl-d23fa78d"><span class="cl-d23fa78b">3.086</span></p></td><td class="cl-d24209dc"><p class="cl-d23fa78d"><span class="cl-d23fa78b">16.000</span></p></td><td class="cl-d24209dd"><p class="cl-d23fa78d"><span class="cl-d23fa78b">51.762</span></p></td><td class="cl-d24209da"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d24209db"><p class="cl-d23fa78c"><span class="cl-d23fa78b">ti(ewood_precip,et)</span></p></td><td class="cl-d24209df"><p class="cl-d23fa78d"><span class="cl-d23fa78b">10.693</span></p></td><td class="cl-d24209dc"><p class="cl-d23fa78d"><span class="cl-d23fa78b">16.000</span></p></td><td class="cl-d24209dd"><p class="cl-d23fa78d"><span class="cl-d23fa78b">390.136</span></p></td><td class="cl-d24209da"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d24209db"><p class="cl-d23fa78c"><span class="cl-d23fa78b">ti(ewood_tmax,wetness)</span></p></td><td class="cl-d24209df"><p class="cl-d23fa78d"><span class="cl-d23fa78b">12.200</span></p></td><td class="cl-d24209dc"><p class="cl-d23fa78d"><span class="cl-d23fa78b">16.000</span></p></td><td class="cl-d24209dd"><p class="cl-d23fa78d"><span class="cl-d23fa78b">398.489</span></p></td><td class="cl-d24209da"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d24209db"><p class="cl-d23fa78c"><span class="cl-d23fa78b">ti(ewood_tmax,et)</span></p></td><td class="cl-d24209df"><p class="cl-d23fa78d"><span class="cl-d23fa78b">7.821</span></p></td><td class="cl-d24209dc"><p class="cl-d23fa78d"><span class="cl-d23fa78b">16.000</span></p></td><td class="cl-d24209dd"><p class="cl-d23fa78d"><span class="cl-d23fa78b">134.327</span></p></td><td class="cl-d24209da"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-d24209e2"><p class="cl-d23fa78c"><span class="cl-d23fa78b">ti(wetness,et)</span></p></td><td class="cl-d2446c0c"><p class="cl-d23fa78d"><span class="cl-d23fa78b">11.091</span></p></td><td class="cl-d24209e1"><p class="cl-d23fa78d"><span class="cl-d23fa78b">16.000</span></p></td><td class="cl-d24209e0"><p class="cl-d23fa78d"><span class="cl-d23fa78b">469.693</span></p></td><td class="cl-d24209e3"><p class="cl-d23fa78d"><span class="cl-d23fa78b">***</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-d246ce59"><p class="cl-d23fa78c"><span class="cl-d23fa78b">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-d246ce59"><p class="cl-d23fa78c"><span class="cl-d23fa78b"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-d246ce59"><p class="cl-d23fa78c"><span class="cl-d23fa78b">Adjusted R-squared: 0.169, Deviance explained 0.683</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-d246ce59"><p class="cl-d23fa78c"><span class="cl-d23fa78b">-REML : 27.805, Scale est: 1.000, N: 312</span></p></td></tr></tfoot></table></div></template>
<div class="flextable-shadow-host" id="9721eea5-b66d-475b-b740-6fff2206bc92"></div>
<script>
var dest = document.getElementById("9721eea5-b66d-475b-b740-6fff2206bc92");
var template = document.getElementById("1e1f21ef-add5-4b6d-a039-3379ea594b8c");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;text-align:center;";
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

```



```r
# we need a hurdle model here. We will fit logistic regression to predict 
# the probabilty of a non-zero flow
# and a scaled-t model with log link to predict the mean of the non-zero data 

m3.gam <- gam(non_zero ~
            s(ewood_precip, bs = "ts") +
            s(ewood_tmax, bs = "ts") +
            s(lagPrecip, bs = "ts") +
            s(doy, bs = "cc") +
            s(wetness, bs = "ts") +
            s(et, bs = "ts") +
            ti(ewood_precip, ewood_tmax) +
            ti(ewood_precip, wetness) +
            ti(ewood_precip, et) +
            ti(ewood_tmax, wetness) +
            ti(ewood_tmax, et) +
            ti(wetness, et),
          data = df_16882,
          select = TRUE,
          family = binomial(),
          method = "REML",
          control = gam.control(nthreads = 2))

m4.gam <- gam(adjusted_flow ~
            #s(ewood_precip, bs = "ts", k = 5) +
            #s(ewood_tmax, bs = "ts", k = 5) +
            s(doy, bs = "cs", k = 6) +
            s(lagPrecip, bs = "ts", k = 5) +
            #s(wetness, bs = "ts", k = 4) +
            #s(et, bs = "ts", k = 4) +
            #te(ewood_precip, ewood_tmax, k = 8),
            te(ewood_precip, wetness, k = 3)+
            te(ewood_precip, et, k = 3) +
            te(ewood_tmax, wetness, k = 3)+
            te(ewood_tmax, et, k = 3) +
            te(wetness, et, k = 3),
          data = subset(df_16882, non_zero == 1),
          select = TRUE,
          family = scat(link = "log"),
          method = "REML",
          control = gam.control(nthreads = 2))
#draw(m4.gam, n = 1000, dist = 2)

flextable::as_flextable(m3.gam) %>%
  flextable::set_caption("Summary of GAMcoefficients at SWQM 16882")
```

```{=html}
<template id="d6a2548d-bb89-40e4-afd5-548657605d83"><style>
.tabwid table{
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 1.275em;
  margin-bottom: 1.275em;
  border-spacing: 0;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-4ae67672{border-collapse:collapse;}.cl-4ad48c8c{font-family:'Arial';font-size:11pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-4ad48c8d{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-4ad48c8e{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4ad48c8f{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4ad48c90{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4ad6ef9a{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad6ef9b{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad6ef9c{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad6ef9d{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad6ef9e{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad6ef9f{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad6efa0{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad6efa1{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad6efa2{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad6efa3{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad6efa4{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad951cc{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad951cd{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad951ce{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad951cf{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad951d0{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad951d1{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad951d2{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad951d3{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad951d4{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad951d5{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4ad951d6{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adbb426{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adbb427{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adbb428{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adbb429{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adbb42a{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adbb42b{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adbb42c{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adbb42d{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adbb42e{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adbb42f{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adbb430{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adcb2fe{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4adcb2ff{width:167.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-4ae67672'>
```
<caption class="Table Caption">(\#tab:gamresults16882)Summary of GAMcoefficients at SWQM 16882</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-4adbb42e"><p class="cl-4ad48c8e"><span class="cl-4ad48c8c">Component</span></p></td><td class="cl-4adcb2ff"><p class="cl-4ad48c8e"><span class="cl-4ad48c8c">Term</span></p></td><td class="cl-4adbb430"><p class="cl-4ad48c8f"><span class="cl-4ad48c8c">Estimate</span></p></td><td class="cl-4adbb42f"><p class="cl-4ad48c8f"><span class="cl-4ad48c8c">Std Error</span></p></td><td class="cl-4adcb2fe"><p class="cl-4ad48c8f"><span class="cl-4ad48c8c">t-value</span></p></td><td class="cl-4adcb2fe"><p class="cl-4ad48c8f"><span class="cl-4ad48c8c">p-value</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-4ad6ef9c"><p class="cl-4ad48c90"><span class="cl-4ad48c8d">A. parametric coefficients</span></p></td><td class="cl-4ad6ef9e"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">(Intercept)</span></p></td><td class="cl-4ad6ef9a"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">-2.527</span></p></td><td class="cl-4ad6ef9d"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.437</span></p></td><td class="cl-4ad6ef9b"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">-5.777</span></p></td><td class="cl-4ad6ef9b"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4ad951d3"><p class="cl-4ad48c90"><span class="cl-4ad48c8c">Component</span></p></td><td class="cl-4ad951d2"><p class="cl-4ad48c8e"><span class="cl-4ad48c8c">Term</span></p></td><td class="cl-4ad951d0"><p class="cl-4ad48c8f"><span class="cl-4ad48c8c">edf</span></p></td><td class="cl-4ad951d4"><p class="cl-4ad48c8f"><span class="cl-4ad48c8c">Ref. df</span></p></td><td class="cl-4ad951d1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8c">F-value</span></p></td><td class="cl-4ad951d1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8c">p-value</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  rowspan="12"class="cl-4adbb428"><p class="cl-4ad48c90"><span class="cl-4ad48c8d">B. smooth terms</span></p></td><td class="cl-4ad951d6"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">s(ewood_precip)</span></p></td><td class="cl-4adbb427"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.927</span></p></td><td class="cl-4adbb426"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">9.000</span></p></td><td class="cl-4ad951d5"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">7.665</span></p></td><td class="cl-4ad951d5"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">**</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4ad6efa3"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">s(ewood_tmax)</span></p></td><td class="cl-4ad6ef9f"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.721</span></p></td><td class="cl-4ad6efa0"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">9.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">2.176</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">.</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4ad6efa3"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">s(lagPrecip)</span></p></td><td class="cl-4ad6ef9f"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.861</span></p></td><td class="cl-4ad6efa0"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">9.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">4.366</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">*</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4ad6efa3"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">s(doy)</span></p></td><td class="cl-4ad6ef9f"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">6.853</span></p></td><td class="cl-4ad6efa0"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">8.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">33.130</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4ad6efa3"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">s(wetness)</span></p></td><td class="cl-4ad6ef9f"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">2.296</span></p></td><td class="cl-4ad6efa0"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">9.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">14.589</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4ad6efa3"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">s(et)</span></p></td><td class="cl-4ad6ef9f"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.324</span></p></td><td class="cl-4ad6efa0"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">9.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.450</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4ad6efa3"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">ti(ewood_precip,ewood_tmax)</span></p></td><td class="cl-4ad6ef9f"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.000</span></p></td><td class="cl-4ad6efa0"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">16.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4ad6efa3"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">ti(ewood_precip,wetness)</span></p></td><td class="cl-4ad6ef9f"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.000</span></p></td><td class="cl-4ad6efa0"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">16.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4ad6efa3"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">ti(ewood_precip,et)</span></p></td><td class="cl-4ad6ef9f"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">1.416</span></p></td><td class="cl-4ad6efa0"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">15.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">5.520</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">*</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4ad6efa3"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">ti(ewood_tmax,wetness)</span></p></td><td class="cl-4ad6ef9f"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.000</span></p></td><td class="cl-4ad6efa0"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">16.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4ad6efa3"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">ti(ewood_tmax,et)</span></p></td><td class="cl-4ad6ef9f"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.000</span></p></td><td class="cl-4ad6efa0"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">16.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">0.000</span></p></td><td class="cl-4ad6efa1"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4ad951ce"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">ti(wetness,et)</span></p></td><td class="cl-4ad6efa4"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">1.769</span></p></td><td class="cl-4ad951cf"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">16.000</span></p></td><td class="cl-4ad951cc"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">4.316</span></p></td><td class="cl-4ad951cc"><p class="cl-4ad48c8f"><span class="cl-4ad48c8d">.</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4adbb42d"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4adbb42d"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4adbb42d"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">Adjusted R-squared: 0.705, Deviance explained 0.672</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4adbb42d"><p class="cl-4ad48c8e"><span class="cl-4ad48c8d">-REML : 90.824, Scale est: 1.000, N: 318</span></p></td></tr></tfoot></table></div></template>
<div class="flextable-shadow-host" id="ced79cdd-40a4-469d-87c4-9dbe278b9304"></div>
<script>
var dest = document.getElementById("ced79cdd-40a4-469d-87c4-9dbe278b9304");
var template = document.getElementById("d6a2548d-bb89-40e4-afd5-548657605d83");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;text-align:center;";
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

```

```r
flextable::as_flextable(m4.gam) %>%
  flextable::set_caption("Summary of GAMcoefficients at SWQM 16882")
```

```{=html}
<template id="d2e92c29-a2c2-4dc6-a428-2029bb983b97"><style>
.tabwid table{
  border-collapse:collapse;
  line-height:1;
  margin-left:auto;
  margin-right:auto;
  border-width: 0;
  display: table;
  margin-top: 1.275em;
  margin-bottom: 1.275em;
  border-spacing: 0;
  border-color: transparent;
}
.tabwid_left table{
  margin-left:0;
}
.tabwid_right table{
  margin-right:0;
}
.tabwid td {
    padding: 0;
}
.tabwid a {
  text-decoration: none;
}
.tabwid thead {
    background-color: transparent;
}
.tabwid tfoot {
    background-color: transparent;
}
.tabwid table tr {
background-color: transparent;
}
</style><div class="tabwid"><style>.cl-4b3abd7c{border-collapse:collapse;}.cl-4b291df6{font-family:'Arial';font-size:11pt;font-weight:bold;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-4b291df7{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-4b291df8{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4b291df9{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4b291dfa{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:5pt;padding-top:5pt;padding-left:5pt;padding-right:5pt;line-height: 1;background-color:transparent;}.cl-4b291dfb{width:148.9pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b291dfc{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b291dfd{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b291dfe{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b291dff{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b291e00{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2b814a{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2b814b{width:148.9pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2b814c{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2b814d{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2b814e{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2b814f{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2b8150{width:148.9pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2b8151{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2b8152{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2b8153{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2b8154{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2de390{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2de391{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2de392{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2de393{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2de394{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2de395{width:148.9pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2de396{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2de397{width:148.9pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2de398{width:144pt;background-color:transparent;vertical-align: top;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2de399{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b2de39a{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b3045ea{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b3045eb{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b3045ec{width:148.9pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b3045ed{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b3045ee{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b304a2c{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b304a2d{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b304a2e{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b304a2f{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b304a30{width:66.4pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b304a31{width:59pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b304a32{width:68.2pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b304a33{width:60.3pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4b304a34{width:148.9pt;background-color:transparent;vertical-align: middle;border-bottom: 1pt solid rgba(102, 102, 102, 1.00);border-top: 1pt solid rgba(102, 102, 102, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-4b3abd7c'>
```
<caption class="Table Caption">(\#tab:gamresults16882)Summary of GAMcoefficients at SWQM 16882</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-4b304a2f"><p class="cl-4b291df8"><span class="cl-4b291df6">Component</span></p></td><td class="cl-4b304a34"><p class="cl-4b291df8"><span class="cl-4b291df6">Term</span></p></td><td class="cl-4b304a30"><p class="cl-4b291df9"><span class="cl-4b291df6">Estimate</span></p></td><td class="cl-4b304a32"><p class="cl-4b291df9"><span class="cl-4b291df6">Std Error</span></p></td><td class="cl-4b304a33"><p class="cl-4b291df9"><span class="cl-4b291df6">t-value</span></p></td><td class="cl-4b304a31"><p class="cl-4b291df9"><span class="cl-4b291df6">p-value</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-4b291e00"><p class="cl-4b291dfa"><span class="cl-4b291df7">A. parametric coefficients</span></p></td><td class="cl-4b291dfb"><p class="cl-4b291df8"><span class="cl-4b291df7">(Intercept)</span></p></td><td class="cl-4b291dff"><p class="cl-4b291df9"><span class="cl-4b291df7">1.370</span></p></td><td class="cl-4b291dfd"><p class="cl-4b291df9"><span class="cl-4b291df7">0.049</span></p></td><td class="cl-4b291dfc"><p class="cl-4b291df9"><span class="cl-4b291df7">27.986</span></p></td><td class="cl-4b291dfe"><p class="cl-4b291df9"><span class="cl-4b291df7">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4b2b814f"><p class="cl-4b291dfa"><span class="cl-4b291df6">Component</span></p></td><td class="cl-4b2b814b"><p class="cl-4b291df8"><span class="cl-4b291df6">Term</span></p></td><td class="cl-4b2b814a"><p class="cl-4b291df9"><span class="cl-4b291df6">edf</span></p></td><td class="cl-4b2b814d"><p class="cl-4b291df9"><span class="cl-4b291df6">Ref. df</span></p></td><td class="cl-4b2b814e"><p class="cl-4b291df9"><span class="cl-4b291df6">F-value</span></p></td><td class="cl-4b2b814c"><p class="cl-4b291df9"><span class="cl-4b291df6">p-value</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  rowspan="7"class="cl-4b2de390"><p class="cl-4b291dfa"><span class="cl-4b291df7">B. smooth terms</span></p></td><td class="cl-4b2b8150"><p class="cl-4b291df8"><span class="cl-4b291df7">s(doy)</span></p></td><td class="cl-4b2b8154"><p class="cl-4b291df9"><span class="cl-4b291df7">4.857</span></p></td><td class="cl-4b2b8152"><p class="cl-4b291df9"><span class="cl-4b291df7">5.000</span></p></td><td class="cl-4b2b8151"><p class="cl-4b291df9"><span class="cl-4b291df7">843.911</span></p></td><td class="cl-4b2b8153"><p class="cl-4b291df9"><span class="cl-4b291df7">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4b2de395"><p class="cl-4b291df8"><span class="cl-4b291df7">s(lagPrecip)</span></p></td><td class="cl-4b2de393"><p class="cl-4b291df9"><span class="cl-4b291df7">2.602</span></p></td><td class="cl-4b2de392"><p class="cl-4b291df9"><span class="cl-4b291df7">4.000</span></p></td><td class="cl-4b2de394"><p class="cl-4b291df9"><span class="cl-4b291df7">66.108</span></p></td><td class="cl-4b2de391"><p class="cl-4b291df9"><span class="cl-4b291df7">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4b2de395"><p class="cl-4b291df8"><span class="cl-4b291df7">te(ewood_precip,wetness)</span></p></td><td class="cl-4b2de393"><p class="cl-4b291df9"><span class="cl-4b291df7">3.252</span></p></td><td class="cl-4b2de392"><p class="cl-4b291df9"><span class="cl-4b291df7">8.000</span></p></td><td class="cl-4b2de394"><p class="cl-4b291df9"><span class="cl-4b291df7">67.099</span></p></td><td class="cl-4b2de391"><p class="cl-4b291df9"><span class="cl-4b291df7">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4b2de395"><p class="cl-4b291df8"><span class="cl-4b291df7">te(ewood_precip,et)</span></p></td><td class="cl-4b2de393"><p class="cl-4b291df9"><span class="cl-4b291df7">3.128</span></p></td><td class="cl-4b2de392"><p class="cl-4b291df9"><span class="cl-4b291df7">6.000</span></p></td><td class="cl-4b2de394"><p class="cl-4b291df9"><span class="cl-4b291df7">456.516</span></p></td><td class="cl-4b2de391"><p class="cl-4b291df9"><span class="cl-4b291df7">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4b2de395"><p class="cl-4b291df8"><span class="cl-4b291df7">te(ewood_tmax,wetness)</span></p></td><td class="cl-4b2de393"><p class="cl-4b291df9"><span class="cl-4b291df7">2.613</span></p></td><td class="cl-4b2de392"><p class="cl-4b291df9"><span class="cl-4b291df7">6.000</span></p></td><td class="cl-4b2de394"><p class="cl-4b291df9"><span class="cl-4b291df7">92.376</span></p></td><td class="cl-4b2de391"><p class="cl-4b291df9"><span class="cl-4b291df7">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4b2de395"><p class="cl-4b291df8"><span class="cl-4b291df7">te(ewood_tmax,et)</span></p></td><td class="cl-4b2de393"><p class="cl-4b291df9"><span class="cl-4b291df7">6.533</span></p></td><td class="cl-4b2de392"><p class="cl-4b291df9"><span class="cl-4b291df7">8.000</span></p></td><td class="cl-4b2de394"><p class="cl-4b291df9"><span class="cl-4b291df7">277.985</span></p></td><td class="cl-4b2de391"><p class="cl-4b291df9"><span class="cl-4b291df7">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-4b2de397"><p class="cl-4b291df8"><span class="cl-4b291df7">te(wetness,et)</span></p></td><td class="cl-4b2de399"><p class="cl-4b291df9"><span class="cl-4b291df7">0.219</span></p></td><td class="cl-4b3045eb"><p class="cl-4b291df9"><span class="cl-4b291df7">4.000</span></p></td><td class="cl-4b2de39a"><p class="cl-4b291df9"><span class="cl-4b291df7">0.334</span></p></td><td class="cl-4b3045ea"><p class="cl-4b291df9"><span class="cl-4b291df7"></span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4b304a2e"><p class="cl-4b291df8"><span class="cl-4b291df7">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4b304a2e"><p class="cl-4b291df8"><span class="cl-4b291df7"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4b304a2e"><p class="cl-4b291df8"><span class="cl-4b291df7">Adjusted R-squared: 0.443, Deviance explained 0.740</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4b304a2e"><p class="cl-4b291df8"><span class="cl-4b291df7">-REML : 223.855, Scale est: 1.000, N: 78</span></p></td></tr></tfoot></table></div></template>
<div class="flextable-shadow-host" id="8a3d87d9-422d-4f86-a56e-df11ee85227e"></div>
<script>
var dest = document.getElementById("8a3d87d9-422d-4f86-a56e-df11ee85227e");
var template = document.getElementById("d2e92c29-a2c2-4dc6-a428-2029bb983b97");
var caption = template.content.querySelector("caption");
if(caption) {
  caption.style.cssText = "display:block;text-align:center;";
  var newcapt = document.createElement("p");
  newcapt.appendChild(caption)
  dest.parentNode.insertBefore(newcapt, dest.previousSibling);
}
var fantome = dest.attachShadow({mode: 'open'});
var templateContent = template.content;
fantome.appendChild(templateContent);
</script>

```



```r
##get model fits and confidence intervals
gam_fam <- family(m1.gam)
ilink <- gam_fam$linkinv

df_16396 %>%
    bind_cols(
    as_tibble(predict(m1.gam, df_16396, se.fit = TRUE)) %>%
      mutate(response = ilink(fit),
             upr_ci = ilink(fit + (2*se.fit)),
             lwr_ci = ilink(fit - (2*se.fit))) %>%
      dplyr::select(response, upr_ci, lwr_ci)
    ) -> df_16396_gam

##get model fits and confidence intervals
gam_fam <- family(m2.gam)
ilink <- gam_fam$linkinv

df_16397 %>%
    bind_cols(
    as_tibble(predict(m2.gam, df_16397, se.fit = TRUE)) %>%
      mutate(response = ilink(fit),
             upr_ci = ilink(fit + (2*se.fit)),
             lwr_ci = ilink(fit - (2*se.fit))) %>%
      dplyr::select(response, upr_ci, lwr_ci)
    ) -> df_16397_gam


##get model fits and confidence intervals
gam_fam <- family(m3.gam)
ilink <- gam_fam$linkinv

df_16882 %>%
    bind_cols(
    as_tibble(predict(m3.gam, df_16882, se.fit = FALSE,
                      type = "response"))) %>%
  mutate(pred_zero = case_when(
    value < 0.5 ~ 0,
    value > 0.5 ~ 1)) %>%
  bind_cols(
    as_tibble(predict(m4.gam, df_16882, se.fit = FALSE,
                      type = "response")) %>%
      dplyr::rename(response = value)) %>%
  mutate(response = case_when(
    pred_zero == 0 ~ 0,
    pred_zero != 0 ~ response)) -> df_16882_gam



df_16396_gam %>%
  ggplot() +
  geom_point(aes(response, adjusted_flow,
                 color = "GAM Estimates against Naturalized Flow"),
             alpha = 0.4) +
  geom_smooth(aes(response, adjusted_flow,
                  color = "GAM Estimates against Naturalized Flow"),
              method = "lm", se = FALSE) +
  geom_abline(aes(linetype = "1:1 Line", intercept = 0, slope = 1)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "GAM Estimated Flow [cfs]", y = "Naturalized Flow [cfs]") +
  guides(color = guide_legend(""), linetype = guide_legend("")) +
  theme_ms() +
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "vertical",
        legend.box.margin = margin(0,0,0,0),
        legend.direction = "vertical",
        legend.position = "none")-> p1

df_16396_gam %>%
  ggplot() +
  geom_line(aes(date, response, color =  "GAM Estimated Flow")) +
  geom_line(aes(date, adjusted_flow, color = "Naturalized Flow SWQM 16396"),
            alpha = 0.4) +
  labs(x = "Date", y = "Mean Daily Flow [cfs]") +
  theme_ms() +
  theme(legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "none") -> p2




df_16397_gam %>%
  ggplot() +
  geom_point(aes(response, adjusted_flow,
                 color = "GAM Estimates against Naturalized Flow"),
             alpha = 0.4) +
  geom_smooth(aes(response, adjusted_flow,
                  color = "GAM Estimates against Naturalized Flow"),
              method = "lm", se = FALSE) +
  geom_abline(aes(linetype = "1:1 Line", intercept = 0, slope = 1)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "GAM Estimated Flow [cfs]", y = "Naturalized Flow [cfs]") +
  guides(color = guide_legend(""), linetype = guide_legend("")) +
  theme_ms() +
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "vertical",
        legend.box.margin = margin(0,0,0,0),
        legend.direction = "vertical",
        legend.position = "none")-> p3

df_16397_gam %>%
  ggplot() +
  geom_line(aes(date, response, color =  "GAM Estimated Flow")) +
  geom_line(aes(date, adjusted_flow, color = "Naturalized Flow SWQM 16397"),
            alpha = 0.4) +
  labs(x = "Date", y = "Mean Daily Flow [cfs]") +
  theme_ms() +
  theme(legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "none") -> p4


df_16882_gam %>%
  ggplot() +
  geom_point(aes(response, adjusted_flow,
                 color = "GAM Estimates against Naturalized Flow"),
             alpha = 0.4) +
  geom_smooth(aes(response, adjusted_flow,
                  color = "GAM Estimates against Naturalized Flow"),
              method = "lm", se = FALSE) +
  geom_abline(aes(linetype = "1:1 Line", intercept = 0, slope = 1)) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "GAM Estimated Flow [cfs]", y = "Naturalized Flow [cfs]") +
  guides(color = guide_legend(""), linetype = guide_legend("")) +
  theme_ms() +
  theme(legend.margin = margin(0,0,0,0),
        legend.box = "vertical",
        legend.box.margin = margin(0,0,0,0),
        legend.direction = "vertical")-> p5

df_16882_gam %>%
  ggplot() +
  geom_line(aes(date, response, color =  "GAM Estimated Flow")) +
  geom_line(aes(date, adjusted_flow, color = "Naturalized Flow"),
            alpha = 0.4) +
  labs(x = "Date", y = "Mean Daily Flow [cfs]") +
  theme_ms() +
  theme(legend.title = element_blank(),
        legend.direction = "vertical") -> p6

(p1+p2)/(p3 +p4)/(p5+p6)
```

<div class="figure">
<img src="document_files/figure-html/gamplots-1.png" alt="GAM estimated streamflows at SWQM 16396 plotted against measured streamflows and over time." width="100%" />
<p class="caption">(\#fig:gamplots)GAM estimated streamflows at SWQM 16396 plotted against measured streamflows and over time.</p>
</div>





```r
## calculate NSE and KGE as goodness of fit metrics
## for each method

##
tibble(model = c("DAR_08065800", "DAR_08109800", "DAR_08110100",
                 "Linear Regression",
                 "GAM"),
       'SWQM-16396' = c(
         hydroGOF::NSE(dar_results_16396$DAR_Q_08065800, as.numeric(dar_results_16396$adjusted_flow)),
         hydroGOF::NSE(dar_results_16396$DAR_Q_08109800, as.numeric(dar_results_16396$adjusted_flow)),
         hydroGOF::NSE(dar_results_16396$DAR_Q_08110100, as.numeric(dar_results_16396$adjusted_flow)),
         hydroGOF::NSE(df_16396_lm$fits, df_16396_lm$adjusted_flow),
         hydroGOF::NSE(as.numeric(df_16396_gam$response), df_16396_gam$adjusted_flow)),
       'SWQM-16397' = c(
         hydroGOF::NSE(dar_results_16397$DAR_Q_08065800, as.numeric(dar_results_16397$adjusted_flow)),
         hydroGOF::NSE(dar_results_16397$DAR_Q_08109800, as.numeric(dar_results_16397$adjusted_flow)),
         hydroGOF::NSE(dar_results_16397$DAR_Q_08110100, as.numeric(dar_results_16397$adjusted_flow)),
         hydroGOF::NSE(df_16397_lm$fits, df_16397_lm$adjusted_flow),
         hydroGOF::NSE(as.numeric(df_16397_gam$response), df_16397_gam$adjusted_flow)),
       'SWQM-16882' = c(
         hydroGOF::NSE(dar_results_16882$DAR_Q_08065800, as.numeric(dar_results_16882$adjusted_flow)),
         hydroGOF::NSE(dar_results_16882$DAR_Q_08109800, as.numeric(dar_results_16882$adjusted_flow)),
         hydroGOF::NSE(dar_results_16882$DAR_Q_08110100, as.numeric(dar_results_16882$adjusted_flow)),
         hydroGOF::NSE(df_16882_lm$fits, df_16882_lm$adjusted_flow),
         hydroGOF::NSE(as.numeric(df_16882_gam$response), df_16882_gam$adjusted_flow)))
```

```
## # A tibble: 5 x 4
##   model             `SWQM-16396` `SWQM-16397` `SWQM-16882`
##   <chr>                    <dbl>        <dbl>        <dbl>
## 1 DAR_08065800           -6.84       -407.         -91.5  
## 2 DAR_08109800            0.0293       -0.462        0.151
## 3 DAR_08110100            0.0996       -7.05        -1.35 
## 4 Linear Regression      -0.603         0.319       -0.237
## 5 GAM                     0.600         0.372       -0.109
```

```r
tibble(model = c("DAR_08065800", "DAR_08109800", "DAR_08110100",
                 "Linear Regression",
                 "GAM"),
       'SWQM-16396' = c(
         hydroGOF::KGE(dar_results_16396$DAR_Q_08065800, as.numeric(dar_results_16396$adjusted_flow)),
         hydroGOF::KGE(dar_results_16396$DAR_Q_08109800, as.numeric(dar_results_16396$adjusted_flow)),
         hydroGOF::KGE(dar_results_16396$DAR_Q_08110100, as.numeric(dar_results_16396$adjusted_flow)),
         hydroGOF::KGE(df_16396_lm$fits, df_16396_lm$adjusted_flow),
         hydroGOF::KGE(as.numeric(df_16396_gam$response), df_16396_gam$adjusted_flow)),
       'SWQM-16397' = c(
         hydroGOF::KGE(dar_results_16397$DAR_Q_08065800, as.numeric(dar_results_16397$adjusted_flow)),
         hydroGOF::KGE(dar_results_16397$DAR_Q_08109800, as.numeric(dar_results_16397$adjusted_flow)),
         hydroGOF::KGE(dar_results_16397$DAR_Q_08110100, as.numeric(dar_results_16397$adjusted_flow)),
         hydroGOF::KGE(df_16397_lm$fits, df_16397_lm$adjusted_flow),
         hydroGOF::KGE(as.numeric(df_16397_gam$response), df_16397_gam$adjusted_flow)),
       'SWQM-16882' = c(
         hydroGOF::KGE(dar_results_16882$DAR_Q_08065800, as.numeric(dar_results_16882$adjusted_flow)),
         hydroGOF::KGE(dar_results_16882$DAR_Q_08109800, as.numeric(dar_results_16882$adjusted_flow)),
         hydroGOF::KGE(dar_results_16882$DAR_Q_08110100, as.numeric(dar_results_16882$adjusted_flow)),
         hydroGOF::KGE(df_16882_lm$fits, df_16882_lm$adjusted_flow),
         hydroGOF::KGE(as.numeric(df_16882_gam$response), df_16882_gam$adjusted_flow)))
```

```
## # A tibble: 5 x 4
##   model             `SWQM-16396` `SWQM-16397` `SWQM-16882`
##   <chr>                    <dbl>        <dbl>        <dbl>
## 1 DAR_08065800           -1.02        -18.4        -9.32  
## 2 DAR_08109800           -0.308         0.201       0.182 
## 3 DAR_08110100           -0.0809       -1.21        0.0916
## 4 Linear Regression       0.277         0.313       0.138 
## 5 GAM                     0.635         0.452       0.402
```

```r
tibble(model = c("DAR_08065800", "DAR_08109800", "DAR_08110100",
                 "Linear Regression",
                 "GAM"),
       'SWQM-16396' = c(
         hydroGOF::pbias(dar_results_16396$DAR_Q_08065800, as.numeric(dar_results_16396$adjusted_flow)),
         hydroGOF::pbias(dar_results_16396$DAR_Q_08109800, as.numeric(dar_results_16396$adjusted_flow)),
         hydroGOF::pbias(dar_results_16396$DAR_Q_08110100, as.numeric(dar_results_16396$adjusted_flow)),
         hydroGOF::pbias(df_16396_lm$fits, df_16396_lm$adjusted_flow),
         hydroGOF::pbias(as.numeric(df_16396_gam$response), df_16396_gam$adjusted_flow)),
       'SWQM-16397' = c(
         hydroGOF::pbias(dar_results_16397$DAR_Q_08065800, as.numeric(dar_results_16397$adjusted_flow)),
         hydroGOF::pbias(dar_results_16397$DAR_Q_08109800, as.numeric(dar_results_16397$adjusted_flow)),
         hydroGOF::pbias(dar_results_16397$DAR_Q_08110100, as.numeric(dar_results_16397$adjusted_flow)),
         hydroGOF::pbias(df_16397_lm$fits, df_16397_lm$adjusted_flow),
         hydroGOF::pbias(as.numeric(df_16397_gam$response), df_16397_gam$adjusted_flow)),
       'SWQM-16882' = c(
         hydroGOF::pbias(dar_results_16882$DAR_Q_08065800, as.numeric(dar_results_16882$adjusted_flow)),
         hydroGOF::pbias(dar_results_16882$DAR_Q_08109800, as.numeric(dar_results_16882$adjusted_flow)),
         hydroGOF::pbias(dar_results_16882$DAR_Q_08110100, as.numeric(dar_results_16882$adjusted_flow)),
         hydroGOF::pbias(df_16882_lm$fits, df_16882_lm$adjusted_flow),
         hydroGOF::pbias(as.numeric(df_16882_gam$response), df_16882_gam$adjusted_flow)))
```

```
## # A tibble: 5 x 4
##   model             `SWQM-16396` `SWQM-16397` `SWQM-16882`
##   <chr>                    <dbl>        <dbl>        <dbl>
## 1 DAR_08065800              58.6        365          581. 
## 2 DAR_08109800             -84.1        -55.9        -33.1
## 3 DAR_08110100             -69.8        -12.1         29.3
## 4 Linear Regression        -22.7         -7.8        -48.6
## 5 GAM                      -15.9         -5.8        -12
```





## References
