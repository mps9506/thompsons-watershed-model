---
title: "Empirical Streamflow Predictions at Thompsons Creek"
date: "2021-04-07"
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
<div class="tabwid"><style>.cl-6423feb6{border-collapse:collapse;}.cl-641a1a68{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-641a1a69{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:italic;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-641a416e{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-641a416f{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-641a8f52{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641a8f53{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641a8f54{width:55pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641a8f55{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641a8f56{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641a8f57{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641a8f58{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641a8f59{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641a8f5a{width:55pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641a8f5b{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641a8f5c{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641ab63a{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641ab63b{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641ab63c{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641ab63d{width:55pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641ab63e{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641ab63f{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641ab640{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641ab641{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641ab642{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641ab643{width:55pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641ab644{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641add2c{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641add2d{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641add2e{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641add2f{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641add30{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641add31{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641add32{width:55pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641add33{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641add34{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641add35{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641add36{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641b050e{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641b050f{width:55pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641b0510{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641b0511{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641b0512{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641b0513{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641b0514{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641b0515{width:55pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-641b0516{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-6423feb6'>
```
<caption class="Table Caption">(\#tab:lrresults16396)Summary of linear regression coefficients at SWQM 16396</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-641b0514"><p class="cl-641a416e"><span class="cl-641a1a68"></span></p></td><td class="cl-641b0512"><p class="cl-641a416f"><span class="cl-641a1a68">Estimate</span></p></td><td class="cl-641b0513"><p class="cl-641a416f"><span class="cl-641a1a68">Standard Error</span></p></td><td class="cl-641b0515"><p class="cl-641a416f"><span class="cl-641a1a68">t value</span></p></td><td class="cl-641b0516"><p class="cl-641a416f"><span class="cl-641a1a68">Pr(&gt;|t|)</span></p></td><td class="cl-641b0511"><p class="cl-641a416e"><span class="cl-641a1a68"></span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-641a8f57"><p class="cl-641a416e"><span class="cl-641a1a68">(Intercept)</span></p></td><td class="cl-641a8f56"><p class="cl-641a416f"><span class="cl-641a1a68">1.83e+00</span></p></td><td class="cl-641a8f55"><p class="cl-641a416f"><span class="cl-641a1a68">4.43e-02</span></p></td><td class="cl-641a8f54"><p class="cl-641a416f"><span class="cl-641a1a68">41.1913</span></p></td><td class="cl-641a8f53"><p class="cl-641a416f"><span class="cl-641a1a68">2.052e-126</span></p></td><td class="cl-641a8f52"><p class="cl-641a416e"><span class="cl-641a1a68">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-641ab63a"><p class="cl-641a416e"><span class="cl-641a1a68">Flow_08065800</span></p></td><td class="cl-641a8f58"><p class="cl-641a416f"><span class="cl-641a1a68">1.79e-04</span></p></td><td class="cl-641a8f59"><p class="cl-641a416f"><span class="cl-641a1a68">9.58e-05</span></p></td><td class="cl-641a8f5a"><p class="cl-641a416f"><span class="cl-641a1a68">1.8656</span></p></td><td class="cl-641a8f5b"><p class="cl-641a416f"><span class="cl-641a1a68">6.305e-02</span></p></td><td class="cl-641a8f5c"><p class="cl-641a416e"><span class="cl-641a1a68">.</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-641ab63a"><p class="cl-641a416e"><span class="cl-641a1a68">Flow_08109800</span></p></td><td class="cl-641a8f58"><p class="cl-641a416f"><span class="cl-641a1a68">1.69e-02</span></p></td><td class="cl-641a8f59"><p class="cl-641a416f"><span class="cl-641a1a68">2.25e-03</span></p></td><td class="cl-641a8f5a"><p class="cl-641a416f"><span class="cl-641a1a68">7.5366</span></p></td><td class="cl-641a8f5b"><p class="cl-641a416f"><span class="cl-641a1a68">5.578e-13</span></p></td><td class="cl-641a8f5c"><p class="cl-641a416e"><span class="cl-641a1a68">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-641ab63a"><p class="cl-641a416e"><span class="cl-641a1a68">Flow_08110100</span></p></td><td class="cl-641a8f58"><p class="cl-641a416f"><span class="cl-641a1a68">5.16e-03</span></p></td><td class="cl-641a8f59"><p class="cl-641a416f"><span class="cl-641a1a68">8.31e-04</span></p></td><td class="cl-641a8f5a"><p class="cl-641a416f"><span class="cl-641a1a68">6.2130</span></p></td><td class="cl-641a8f5b"><p class="cl-641a416f"><span class="cl-641a1a68">1.708e-09</span></p></td><td class="cl-641a8f5c"><p class="cl-641a416e"><span class="cl-641a1a68">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-641ab640"><p class="cl-641a416e"><span class="cl-641a1a68">lag_Flow_08065800</span></p></td><td class="cl-641ab63f"><p class="cl-641a416f"><span class="cl-641a1a68">7.85e-06</span></p></td><td class="cl-641ab63e"><p class="cl-641a416f"><span class="cl-641a1a68">9.50e-05</span></p></td><td class="cl-641ab63d"><p class="cl-641a416f"><span class="cl-641a1a68">0.0827</span></p></td><td class="cl-641ab63c"><p class="cl-641a416f"><span class="cl-641a1a68">9.342e-01</span></p></td><td class="cl-641ab63b"><p class="cl-641a416e"><span class="cl-641a1a68"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-641ab640"><p class="cl-641a416e"><span class="cl-641a1a68">lag_Flow_08109800</span></p></td><td class="cl-641ab63f"><p class="cl-641a416f"><span class="cl-641a1a68">-4.97e-03</span></p></td><td class="cl-641ab63e"><p class="cl-641a416f"><span class="cl-641a1a68">2.23e-03</span></p></td><td class="cl-641ab63d"><p class="cl-641a416f"><span class="cl-641a1a68">-2.2311</span></p></td><td class="cl-641ab63c"><p class="cl-641a416f"><span class="cl-641a1a68">2.640e-02</span></p></td><td class="cl-641ab63b"><p class="cl-641a416e"><span class="cl-641a1a68">*</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-641ab644"><p class="cl-641a416e"><span class="cl-641a1a68">lag_Flow_08110100</span></p></td><td class="cl-641add2d"><p class="cl-641a416f"><span class="cl-641a1a68">-1.76e-03</span></p></td><td class="cl-641add2c"><p class="cl-641a416f"><span class="cl-641a1a68">8.41e-04</span></p></td><td class="cl-641ab643"><p class="cl-641a416f"><span class="cl-641a1a68">-2.0984</span></p></td><td class="cl-641ab642"><p class="cl-641a416f"><span class="cl-641a1a68">3.669e-02</span></p></td><td class="cl-641ab641"><p class="cl-641a416e"><span class="cl-641a1a68">*</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-641add33"><p class="cl-641a416f"><span class="cl-641a1a69">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-641add36"><p class="cl-641a416e"><span class="cl-641a1a68"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-641add36"><p class="cl-641a416e"><span class="cl-641a1a68">Residual standard error: 0.7083 on 304 degrees of freedom</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-641add36"><p class="cl-641a416e"><span class="cl-641a1a68">Multiple R-squared: 0.4088, Adjusted R-squared: 0.3972</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-641add36"><p class="cl-641a416e"><span class="cl-641a1a68">F-statistic: 35.04 on 304 and 6 DF, p-value: 0.0000</span></p></td></tr></tfoot></table></div>
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
<div class="tabwid"><style>.cl-6470d9c0{border-collapse:collapse;}.cl-6466cea8{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-6466cea9{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:italic;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-6466ceaa{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-6466ceab{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-64671ca0{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64671ca1{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64671ca2{width:49pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64671ca3{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64671ca4{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64671ca5{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64671ca6{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64671ca7{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64671ca8{width:49pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64671ca9{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64671caa{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-6467439c{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-6467439d{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-6467439e{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-6467439f{width:49pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-646743a0{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-646743a1{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-646743a2{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-646743a3{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-646743a4{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-646743a5{width:49pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-646743a6{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64676a84{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64676a85{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64676a86{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64676a87{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64676a88{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64676a89{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64676a8a{width:49pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64676a8b{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64676a8c{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64676a8d{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64676a8e{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-6467916c{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-6467916d{width:49pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-6467916e{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-6467916f{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64679170{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64679171{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64679172{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64679173{width:49pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64679174{width:71pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-6470d9c0'>
```
<caption class="Table Caption">(\#tab:lrresults16397)Summary of linear regression coefficients at SWQM 16397</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-64679172"><p class="cl-6466ceaa"><span class="cl-6466cea8"></span></p></td><td class="cl-64679170"><p class="cl-6466ceab"><span class="cl-6466cea8">Estimate</span></p></td><td class="cl-64679171"><p class="cl-6466ceab"><span class="cl-6466cea8">Standard Error</span></p></td><td class="cl-64679173"><p class="cl-6466ceab"><span class="cl-6466cea8">t value</span></p></td><td class="cl-64679174"><p class="cl-6466ceab"><span class="cl-6466cea8">Pr(&gt;|t|)</span></p></td><td class="cl-6467916f"><p class="cl-6466ceaa"><span class="cl-6466cea8"></span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-64671ca5"><p class="cl-6466ceaa"><span class="cl-6466cea8">(Intercept)</span></p></td><td class="cl-64671ca4"><p class="cl-6466ceab"><span class="cl-6466cea8">8.46e-01</span></p></td><td class="cl-64671ca3"><p class="cl-6466ceab"><span class="cl-6466cea8">1.48e-02</span></p></td><td class="cl-64671ca2"><p class="cl-6466ceab"><span class="cl-6466cea8">57.074</span></p></td><td class="cl-64671ca1"><p class="cl-6466ceab"><span class="cl-6466cea8">1.691e-164</span></p></td><td class="cl-64671ca0"><p class="cl-6466ceaa"><span class="cl-6466cea8">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-6467439c"><p class="cl-6466ceaa"><span class="cl-6466cea8">Flow_08065800</span></p></td><td class="cl-64671ca6"><p class="cl-6466ceab"><span class="cl-6466cea8">2.23e-05</span></p></td><td class="cl-64671ca7"><p class="cl-6466ceab"><span class="cl-6466cea8">3.20e-05</span></p></td><td class="cl-64671ca8"><p class="cl-6466ceab"><span class="cl-6466cea8">0.695</span></p></td><td class="cl-64671ca9"><p class="cl-6466ceab"><span class="cl-6466cea8">4.877e-01</span></p></td><td class="cl-64671caa"><p class="cl-6466ceaa"><span class="cl-6466cea8"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-6467439c"><p class="cl-6466ceaa"><span class="cl-6466cea8">Flow_08109800</span></p></td><td class="cl-64671ca6"><p class="cl-6466ceab"><span class="cl-6466cea8">5.64e-03</span></p></td><td class="cl-64671ca7"><p class="cl-6466ceab"><span class="cl-6466cea8">7.51e-04</span></p></td><td class="cl-64671ca8"><p class="cl-6466ceab"><span class="cl-6466cea8">7.509</span></p></td><td class="cl-64671ca9"><p class="cl-6466ceab"><span class="cl-6466cea8">6.650e-13</span></p></td><td class="cl-64671caa"><p class="cl-6466ceaa"><span class="cl-6466cea8">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-6467439c"><p class="cl-6466ceaa"><span class="cl-6466cea8">Flow_08110100</span></p></td><td class="cl-64671ca6"><p class="cl-6466ceab"><span class="cl-6466cea8">2.30e-03</span></p></td><td class="cl-64671ca7"><p class="cl-6466ceab"><span class="cl-6466cea8">2.78e-04</span></p></td><td class="cl-64671ca8"><p class="cl-6466ceab"><span class="cl-6466cea8">8.296</span></p></td><td class="cl-64671ca9"><p class="cl-6466ceab"><span class="cl-6466cea8">3.547e-15</span></p></td><td class="cl-64671caa"><p class="cl-6466ceaa"><span class="cl-6466cea8">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-646743a2"><p class="cl-6466ceaa"><span class="cl-6466cea8">lag_Flow_08065800</span></p></td><td class="cl-646743a1"><p class="cl-6466ceab"><span class="cl-6466cea8">5.52e-05</span></p></td><td class="cl-646743a0"><p class="cl-6466ceab"><span class="cl-6466cea8">3.17e-05</span></p></td><td class="cl-6467439f"><p class="cl-6466ceab"><span class="cl-6466cea8">1.739</span></p></td><td class="cl-6467439e"><p class="cl-6466ceab"><span class="cl-6466cea8">8.304e-02</span></p></td><td class="cl-6467439d"><p class="cl-6466ceaa"><span class="cl-6466cea8">.</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-646743a2"><p class="cl-6466ceaa"><span class="cl-6466cea8">lag_Flow_08109800</span></p></td><td class="cl-646743a1"><p class="cl-6466ceab"><span class="cl-6466cea8">-2.96e-03</span></p></td><td class="cl-646743a0"><p class="cl-6466ceab"><span class="cl-6466cea8">7.45e-04</span></p></td><td class="cl-6467439f"><p class="cl-6466ceab"><span class="cl-6466cea8">-3.970</span></p></td><td class="cl-6467439e"><p class="cl-6466ceab"><span class="cl-6466cea8">8.990e-05</span></p></td><td class="cl-6467439d"><p class="cl-6466ceaa"><span class="cl-6466cea8">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-646743a6"><p class="cl-6466ceaa"><span class="cl-6466cea8">lag_Flow_08110100</span></p></td><td class="cl-64676a85"><p class="cl-6466ceab"><span class="cl-6466cea8">-8.21e-04</span></p></td><td class="cl-64676a84"><p class="cl-6466ceab"><span class="cl-6466cea8">2.81e-04</span></p></td><td class="cl-646743a5"><p class="cl-6466ceab"><span class="cl-6466cea8">-2.920</span></p></td><td class="cl-646743a4"><p class="cl-6466ceab"><span class="cl-6466cea8">3.760e-03</span></p></td><td class="cl-646743a3"><p class="cl-6466ceaa"><span class="cl-6466cea8">**</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-64676a8b"><p class="cl-6466ceab"><span class="cl-6466cea9">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-64676a8e"><p class="cl-6466ceaa"><span class="cl-6466cea8"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-64676a8e"><p class="cl-6466ceaa"><span class="cl-6466cea8">Residual standard error: 0.2368 on 304 degrees of freedom</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-64676a8e"><p class="cl-6466ceaa"><span class="cl-6466cea8">Multiple R-squared: 0.4337, Adjusted R-squared: 0.4226</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-64676a8e"><p class="cl-6466ceaa"><span class="cl-6466cea8">F-statistic: 38.81 on 304 and 6 DF, p-value: 0.0000</span></p></td></tr></tfoot></table></div>
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
<div class="tabwid"><style>.cl-64bb584c{border-collapse:collapse;}.cl-64b3208c{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-64b3208d{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:italic;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-64b34792{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-64b34793{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-64b36e7a{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b36e7b{width:65pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b36e7c{width:53pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b36e7d{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b36e7e{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b36e7f{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b36e80{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b36e81{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b36e82{width:53pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b36e83{width:65pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b36e84{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3959e{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3959f{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b395a0{width:65pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b395a1{width:53pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b395a2{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b395a3{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b395a4{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b395a5{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b395a6{width:65pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b395a7{width:53pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b395a8{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3bc68{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3bc69{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3bc6a{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3bc6b{width:65pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3bc6c{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3bc6d{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3bc6e{width:53pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3bc6f{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3bc70{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3bc71{width:65pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3bc72{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3e36e{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3e36f{width:53pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3e370{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3e371{width:28pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3e372{width:62pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3e373{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3e374{width:114pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3e375{width:53pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-64b3e376{width:65pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-64bb584c'>
```
<caption class="Table Caption">(\#tab:lrresults16882)Summary of linear regression coefficients at SWQM 16882</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-64b3e374"><p class="cl-64b34792"><span class="cl-64b3208c"></span></p></td><td class="cl-64b3e372"><p class="cl-64b34793"><span class="cl-64b3208c">Estimate</span></p></td><td class="cl-64b3e373"><p class="cl-64b34793"><span class="cl-64b3208c">Standard Error</span></p></td><td class="cl-64b3e375"><p class="cl-64b34793"><span class="cl-64b3208c">t value</span></p></td><td class="cl-64b3e376"><p class="cl-64b34793"><span class="cl-64b3208c">Pr(&gt;|t|)</span></p></td><td class="cl-64b3e371"><p class="cl-64b34792"><span class="cl-64b3208c"></span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-64b36e7f"><p class="cl-64b34792"><span class="cl-64b3208c">(Intercept)</span></p></td><td class="cl-64b36e7e"><p class="cl-64b34793"><span class="cl-64b3208c">3.12e-01</span></p></td><td class="cl-64b36e7d"><p class="cl-64b34793"><span class="cl-64b3208c">0.049030</span></p></td><td class="cl-64b36e7c"><p class="cl-64b34793"><span class="cl-64b3208c">6.3536</span></p></td><td class="cl-64b36e7b"><p class="cl-64b34793"><span class="cl-64b3208c">7.497e-10</span></p></td><td class="cl-64b36e7a"><p class="cl-64b34792"><span class="cl-64b3208c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-64b3959e"><p class="cl-64b34792"><span class="cl-64b3208c">Flow_08065800</span></p></td><td class="cl-64b36e80"><p class="cl-64b34793"><span class="cl-64b3208c">-9.30e-06</span></p></td><td class="cl-64b36e81"><p class="cl-64b34793"><span class="cl-64b3208c">0.000107</span></p></td><td class="cl-64b36e82"><p class="cl-64b34793"><span class="cl-64b3208c">-0.0868</span></p></td><td class="cl-64b36e83"><p class="cl-64b34793"><span class="cl-64b3208c">9.309e-01</span></p></td><td class="cl-64b36e84"><p class="cl-64b34792"><span class="cl-64b3208c"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-64b3959e"><p class="cl-64b34792"><span class="cl-64b3208c">Flow_08109800</span></p></td><td class="cl-64b36e80"><p class="cl-64b34793"><span class="cl-64b3208c">1.81e-02</span></p></td><td class="cl-64b36e81"><p class="cl-64b34793"><span class="cl-64b3208c">0.002513</span></p></td><td class="cl-64b36e82"><p class="cl-64b34793"><span class="cl-64b3208c">7.1853</span></p></td><td class="cl-64b36e83"><p class="cl-64b34793"><span class="cl-64b3208c">5.022e-12</span></p></td><td class="cl-64b36e84"><p class="cl-64b34792"><span class="cl-64b3208c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-64b3959e"><p class="cl-64b34792"><span class="cl-64b3208c">Flow_08110100</span></p></td><td class="cl-64b36e80"><p class="cl-64b34793"><span class="cl-64b3208c">3.65e-03</span></p></td><td class="cl-64b36e81"><p class="cl-64b34793"><span class="cl-64b3208c">0.000929</span></p></td><td class="cl-64b36e82"><p class="cl-64b34793"><span class="cl-64b3208c">3.9333</span></p></td><td class="cl-64b36e83"><p class="cl-64b34793"><span class="cl-64b3208c">1.034e-04</span></p></td><td class="cl-64b36e84"><p class="cl-64b34792"><span class="cl-64b3208c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-64b395a4"><p class="cl-64b34792"><span class="cl-64b3208c">lag_Flow_08065800</span></p></td><td class="cl-64b395a3"><p class="cl-64b34793"><span class="cl-64b3208c">6.80e-05</span></p></td><td class="cl-64b395a2"><p class="cl-64b34793"><span class="cl-64b3208c">0.000106</span></p></td><td class="cl-64b395a1"><p class="cl-64b34793"><span class="cl-64b3208c">0.6405</span></p></td><td class="cl-64b395a0"><p class="cl-64b34793"><span class="cl-64b3208c">5.223e-01</span></p></td><td class="cl-64b3959f"><p class="cl-64b34792"><span class="cl-64b3208c"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-64b395a4"><p class="cl-64b34792"><span class="cl-64b3208c">lag_Flow_08109800</span></p></td><td class="cl-64b395a3"><p class="cl-64b34793"><span class="cl-64b3208c">-8.79e-03</span></p></td><td class="cl-64b395a2"><p class="cl-64b34793"><span class="cl-64b3208c">0.002494</span></p></td><td class="cl-64b395a1"><p class="cl-64b34793"><span class="cl-64b3208c">-3.5244</span></p></td><td class="cl-64b395a0"><p class="cl-64b34793"><span class="cl-64b3208c">4.885e-04</span></p></td><td class="cl-64b3959f"><p class="cl-64b34792"><span class="cl-64b3208c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-64b395a8"><p class="cl-64b34792"><span class="cl-64b3208c">lag_Flow_08110100</span></p></td><td class="cl-64b3bc69"><p class="cl-64b34793"><span class="cl-64b3208c">-2.54e-03</span></p></td><td class="cl-64b3bc68"><p class="cl-64b34793"><span class="cl-64b3208c">0.000940</span></p></td><td class="cl-64b395a7"><p class="cl-64b34793"><span class="cl-64b3208c">-2.7016</span></p></td><td class="cl-64b395a6"><p class="cl-64b34793"><span class="cl-64b3208c">7.280e-03</span></p></td><td class="cl-64b395a5"><p class="cl-64b34792"><span class="cl-64b3208c">**</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-64b3bc6f"><p class="cl-64b34793"><span class="cl-64b3208d">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-64b3bc72"><p class="cl-64b34792"><span class="cl-64b3208c"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-64b3bc72"><p class="cl-64b34792"><span class="cl-64b3208c">Residual standard error: 0.7923 on 310 degrees of freedom</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-64b3bc72"><p class="cl-64b34792"><span class="cl-64b3208c">Multiple R-squared: 0.2595, Adjusted R-squared: 0.2452</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-64b3bc72"><p class="cl-64b34792"><span class="cl-64b3208c">F-statistic: 18.11 on 310 and 6 DF, p-value: 0.0000</span></p></td></tr></tfoot></table></div>
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
<div class="tabwid"><style>.cl-55cd44b4{border-collapse:collapse;}.cl-55c14f92{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-55c14f93{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:italic;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-55c29ce4{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-55c29ce5{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-55c29faa{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c29fab{width:50pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c29fac{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c29fad{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c29fae{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c29faf{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c29fb0{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c29fb1{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c29fb2{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c29fb3{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c29fb4{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c429f6{width:50pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c429f7{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c42aa0{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c42aa1{width:50pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c42aa2{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c42aa3{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c42aa4{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c42aa5{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c42aa6{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c42aa7{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c42aa8{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c42aa9{width:50pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c42aaa{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c4ee22{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c4ee23{width:50pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c4ee24{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c4ee25{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c4ee26{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-55c4ee27{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-55cd44b4'>
```
<caption class="Table Caption">(\#tab:gamresults16396)Summary of GAM coefficients at SWQM 16396</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-55c4ee25"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c4ee22"><p class="cl-55c29ce4"><span class="cl-55c14f92">Estimate</span></p></td><td class="cl-55c4ee27"><p class="cl-55c29ce4"><span class="cl-55c14f92">Standard Error</span></p></td><td class="cl-55c4ee23"><p class="cl-55c29ce5"><span class="cl-55c14f92">z value</span></p></td><td class="cl-55c4ee24"><p class="cl-55c29ce5"><span class="cl-55c14f92">Pr(&gt;|z|)</span></p></td><td class="cl-55c4ee26"><p class="cl-55c29ce4"><span class="cl-55c14f92">Signif.</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-55c29fac"><p class="cl-55c29ce4"><span class="cl-55c14f92">s(ewood_precip)</span></p></td><td class="cl-55c29fae"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29faa"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fab"><p class="cl-55c29ce5"><span class="cl-55c14f92">8.17</span></p></td><td class="cl-55c29faf"><p class="cl-55c29ce5"><span class="cl-55c14f92">0.000e+00</span></p></td><td class="cl-55c29fad"><p class="cl-55c29ce4"><span class="cl-55c14f92">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-55c29fac"><p class="cl-55c29ce4"><span class="cl-55c14f92">s(ewood_tmax)</span></p></td><td class="cl-55c29fae"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29faa"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fab"><p class="cl-55c29ce5"><span class="cl-55c14f92">9.40</span></p></td><td class="cl-55c29faf"><p class="cl-55c29ce5"><span class="cl-55c14f92">6.247e-06</span></p></td><td class="cl-55c29fad"><p class="cl-55c29ce4"><span class="cl-55c14f92">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-55c29fac"><p class="cl-55c29ce4"><span class="cl-55c14f92">s(lagPrecip)</span></p></td><td class="cl-55c29fae"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29faa"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fab"><p class="cl-55c29ce5"><span class="cl-55c14f92">342.94</span></p></td><td class="cl-55c29faf"><p class="cl-55c29ce5"><span class="cl-55c14f92">0.000e+00</span></p></td><td class="cl-55c29fad"><p class="cl-55c29ce4"><span class="cl-55c14f92">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-55c29fac"><p class="cl-55c29ce4"><span class="cl-55c14f92">s(doy)</span></p></td><td class="cl-55c29fae"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29faa"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fab"><p class="cl-55c29ce5"><span class="cl-55c14f92">333.11</span></p></td><td class="cl-55c29faf"><p class="cl-55c29ce5"><span class="cl-55c14f92">0.000e+00</span></p></td><td class="cl-55c29fad"><p class="cl-55c29ce4"><span class="cl-55c14f92">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-55c29fac"><p class="cl-55c29ce4"><span class="cl-55c14f92">s(wetness)</span></p></td><td class="cl-55c29fae"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29faa"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fab"><p class="cl-55c29ce5"><span class="cl-55c14f92">50.00</span></p></td><td class="cl-55c29faf"><p class="cl-55c29ce5"><span class="cl-55c14f92">0.000e+00</span></p></td><td class="cl-55c29fad"><p class="cl-55c29ce4"><span class="cl-55c14f92">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-55c29fac"><p class="cl-55c29ce4"><span class="cl-55c14f92">s(et)</span></p></td><td class="cl-55c29fae"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29faa"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fab"><p class="cl-55c29ce5"><span class="cl-55c14f92">1.98</span></p></td><td class="cl-55c29faf"><p class="cl-55c29ce5"><span class="cl-55c14f92">2.336e-02</span></p></td><td class="cl-55c29fad"><p class="cl-55c29ce4"><span class="cl-55c14f92">*</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-55c29fac"><p class="cl-55c29ce4"><span class="cl-55c14f92">ti(ewood_precip,lagPrecip)</span></p></td><td class="cl-55c29fae"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29faa"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fab"><p class="cl-55c29ce5"><span class="cl-55c14f92">19.35</span></p></td><td class="cl-55c29faf"><p class="cl-55c29ce5"><span class="cl-55c14f92">7.282e-07</span></p></td><td class="cl-55c29fad"><p class="cl-55c29ce4"><span class="cl-55c14f92">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-55c29fac"><p class="cl-55c29ce4"><span class="cl-55c14f92">ti(ewood_precip,ewood_tmax)</span></p></td><td class="cl-55c29fae"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29faa"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fab"><p class="cl-55c29ce5"><span class="cl-55c14f92">130.02</span></p></td><td class="cl-55c29faf"><p class="cl-55c29ce5"><span class="cl-55c14f92">0.000e+00</span></p></td><td class="cl-55c29fad"><p class="cl-55c29ce4"><span class="cl-55c14f92">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-55c29fac"><p class="cl-55c29ce4"><span class="cl-55c14f92">ti(ewood_precip,wetness)</span></p></td><td class="cl-55c29fae"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29faa"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fab"><p class="cl-55c29ce5"><span class="cl-55c14f92">139.59</span></p></td><td class="cl-55c29faf"><p class="cl-55c29ce5"><span class="cl-55c14f92">0.000e+00</span></p></td><td class="cl-55c29fad"><p class="cl-55c29ce4"><span class="cl-55c14f92">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-55c29fac"><p class="cl-55c29ce4"><span class="cl-55c14f92">ti(ewood_precip,et)</span></p></td><td class="cl-55c29fae"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29faa"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fab"><p class="cl-55c29ce5"><span class="cl-55c14f92">779.63</span></p></td><td class="cl-55c29faf"><p class="cl-55c29ce5"><span class="cl-55c14f92">0.000e+00</span></p></td><td class="cl-55c29fad"><p class="cl-55c29ce4"><span class="cl-55c14f92">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-55c29fac"><p class="cl-55c29ce4"><span class="cl-55c14f92">ti(ewood_tmax,wetness)</span></p></td><td class="cl-55c29fae"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29faa"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fab"><p class="cl-55c29ce5"><span class="cl-55c14f92">32.16</span></p></td><td class="cl-55c29faf"><p class="cl-55c29ce5"><span class="cl-55c14f92">0.000e+00</span></p></td><td class="cl-55c29fad"><p class="cl-55c29ce4"><span class="cl-55c14f92">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-55c29fac"><p class="cl-55c29ce4"><span class="cl-55c14f92">ti(ewood_tmax,et)</span></p></td><td class="cl-55c29fae"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29faa"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fab"><p class="cl-55c29ce5"><span class="cl-55c14f92">40.48</span></p></td><td class="cl-55c29faf"><p class="cl-55c29ce5"><span class="cl-55c14f92">0.000e+00</span></p></td><td class="cl-55c29fad"><p class="cl-55c29ce4"><span class="cl-55c14f92">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-55c29fb3"><p class="cl-55c29ce4"><span class="cl-55c14f92">ti(wetness,et)</span></p></td><td class="cl-55c29fb4"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c29fb0"><p class="cl-55c29ce4"><span class="cl-55c14f92"></span></p></td><td class="cl-55c429f6"><p class="cl-55c29ce5"><span class="cl-55c14f92">451.08</span></p></td><td class="cl-55c29fb2"><p class="cl-55c29ce5"><span class="cl-55c14f92">0.000e+00</span></p></td><td class="cl-55c29fb1"><p class="cl-55c29ce4"><span class="cl-55c14f92">***</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-55c42aa4"><p class="cl-55c29ce5"><span class="cl-55c14f93">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-55c42aaa"><p class="cl-55c29ce4"><span class="cl-55c14f92"> </span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-55c42aaa"><p class="cl-55c29ce4"><span class="cl-55c14f92">(Dispersion parameter for Scaled t(3.001,1.175) family taken to be 1)</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-55c42aaa"><p class="cl-55c29ce4"><span class="cl-55c14f92">Null deviance: NULL on NULL degrees of freedom</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-55c42aaa"><p class="cl-55c29ce4"><span class="cl-55c14f92">Residual deviance: NULL on NULL degrees of freedom</span></p></td></tr></tfoot></table></div>
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
<div class="tabwid"><style>.cl-b6b76892{border-collapse:collapse;}.cl-b6ab320c{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-b6ab320d{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:italic;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-b6acb960{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-b6acb961{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-b6acb962{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6acb963{width:50pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6acb964{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6acb965{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6acb966{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6acb967{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6acb968{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6acb969{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6acb96a{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3dbc{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3dbd{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3dbe{width:50pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3dbf{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3dc0{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3dc1{width:50pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3dc2{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3dc3{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3dc4{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3e8e{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3e8f{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3e90{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3e91{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3e92{width:50pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3e93{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3e94{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3e95{width:50pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3e96{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3e97{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6ae3e98{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-b6afc7ae{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-b6b76892'>
```
<caption class="Table Caption">(\#tab:gamresults16397)Summary of GAM coefficients at SWQM 16397</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-b6ae3e97"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6ae3e94"><p class="cl-b6acb960"><span class="cl-b6ab320c">Estimate</span></p></td><td class="cl-b6afc7ae"><p class="cl-b6acb960"><span class="cl-b6ab320c">Standard Error</span></p></td><td class="cl-b6ae3e95"><p class="cl-b6acb961"><span class="cl-b6ab320c">z value</span></p></td><td class="cl-b6ae3e96"><p class="cl-b6acb961"><span class="cl-b6ab320c">Pr(&gt;|z|)</span></p></td><td class="cl-b6ae3e98"><p class="cl-b6acb960"><span class="cl-b6ab320c">Signif.</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-b6acb964"><p class="cl-b6acb960"><span class="cl-b6ab320c">s(ewood_precip)</span></p></td><td class="cl-b6acb966"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb962"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb963"><p class="cl-b6acb961"><span class="cl-b6ab320c">11.93</span></p></td><td class="cl-b6acb967"><p class="cl-b6acb961"><span class="cl-b6ab320c">1.229e-06</span></p></td><td class="cl-b6acb965"><p class="cl-b6acb960"><span class="cl-b6ab320c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6acb964"><p class="cl-b6acb960"><span class="cl-b6ab320c">s(ewood_tmax)</span></p></td><td class="cl-b6acb966"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb962"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb963"><p class="cl-b6acb961"><span class="cl-b6ab320c">1.56</span></p></td><td class="cl-b6acb967"><p class="cl-b6acb961"><span class="cl-b6ab320c">4.402e-02</span></p></td><td class="cl-b6acb965"><p class="cl-b6acb960"><span class="cl-b6ab320c">*</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6acb964"><p class="cl-b6acb960"><span class="cl-b6ab320c">s(lagPrecip)</span></p></td><td class="cl-b6acb966"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb962"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb963"><p class="cl-b6acb961"><span class="cl-b6ab320c">54.56</span></p></td><td class="cl-b6acb967"><p class="cl-b6acb961"><span class="cl-b6ab320c">0.000e+00</span></p></td><td class="cl-b6acb965"><p class="cl-b6acb960"><span class="cl-b6ab320c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6acb964"><p class="cl-b6acb960"><span class="cl-b6ab320c">s(doy)</span></p></td><td class="cl-b6acb966"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb962"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb963"><p class="cl-b6acb961"><span class="cl-b6ab320c">285.47</span></p></td><td class="cl-b6acb967"><p class="cl-b6acb961"><span class="cl-b6ab320c">0.000e+00</span></p></td><td class="cl-b6acb965"><p class="cl-b6acb960"><span class="cl-b6ab320c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6acb964"><p class="cl-b6acb960"><span class="cl-b6ab320c">s(wetness)</span></p></td><td class="cl-b6acb966"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb962"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb963"><p class="cl-b6acb961"><span class="cl-b6ab320c">32.63</span></p></td><td class="cl-b6acb967"><p class="cl-b6acb961"><span class="cl-b6ab320c">0.000e+00</span></p></td><td class="cl-b6acb965"><p class="cl-b6acb960"><span class="cl-b6ab320c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6acb964"><p class="cl-b6acb960"><span class="cl-b6ab320c">s(et)</span></p></td><td class="cl-b6acb966"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb962"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb963"><p class="cl-b6acb961"><span class="cl-b6ab320c">26.11</span></p></td><td class="cl-b6acb967"><p class="cl-b6acb961"><span class="cl-b6ab320c">0.000e+00</span></p></td><td class="cl-b6acb965"><p class="cl-b6acb960"><span class="cl-b6ab320c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6acb964"><p class="cl-b6acb960"><span class="cl-b6ab320c">ti(ewood_precip,lagPrecip)</span></p></td><td class="cl-b6acb966"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb962"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb963"><p class="cl-b6acb961"><span class="cl-b6ab320c">209.86</span></p></td><td class="cl-b6acb967"><p class="cl-b6acb961"><span class="cl-b6ab320c">0.000e+00</span></p></td><td class="cl-b6acb965"><p class="cl-b6acb960"><span class="cl-b6ab320c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6acb964"><p class="cl-b6acb960"><span class="cl-b6ab320c">ti(ewood_precip,ewood_tmax)</span></p></td><td class="cl-b6acb966"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb962"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb963"><p class="cl-b6acb961"><span class="cl-b6ab320c">161.10</span></p></td><td class="cl-b6acb967"><p class="cl-b6acb961"><span class="cl-b6ab320c">0.000e+00</span></p></td><td class="cl-b6acb965"><p class="cl-b6acb960"><span class="cl-b6ab320c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6acb964"><p class="cl-b6acb960"><span class="cl-b6ab320c">ti(ewood_precip,wetness)</span></p></td><td class="cl-b6acb966"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb962"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb963"><p class="cl-b6acb961"><span class="cl-b6ab320c">53.93</span></p></td><td class="cl-b6acb967"><p class="cl-b6acb961"><span class="cl-b6ab320c">0.000e+00</span></p></td><td class="cl-b6acb965"><p class="cl-b6acb960"><span class="cl-b6ab320c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6acb964"><p class="cl-b6acb960"><span class="cl-b6ab320c">ti(ewood_precip,et)</span></p></td><td class="cl-b6acb966"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb962"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb963"><p class="cl-b6acb961"><span class="cl-b6ab320c">378.14</span></p></td><td class="cl-b6acb967"><p class="cl-b6acb961"><span class="cl-b6ab320c">0.000e+00</span></p></td><td class="cl-b6acb965"><p class="cl-b6acb960"><span class="cl-b6ab320c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6acb964"><p class="cl-b6acb960"><span class="cl-b6ab320c">ti(ewood_tmax,wetness)</span></p></td><td class="cl-b6acb966"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb962"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb963"><p class="cl-b6acb961"><span class="cl-b6ab320c">375.35</span></p></td><td class="cl-b6acb967"><p class="cl-b6acb961"><span class="cl-b6ab320c">0.000e+00</span></p></td><td class="cl-b6acb965"><p class="cl-b6acb960"><span class="cl-b6ab320c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6acb964"><p class="cl-b6acb960"><span class="cl-b6ab320c">ti(ewood_tmax,et)</span></p></td><td class="cl-b6acb966"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb962"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb963"><p class="cl-b6acb961"><span class="cl-b6ab320c">128.74</span></p></td><td class="cl-b6acb967"><p class="cl-b6acb961"><span class="cl-b6ab320c">0.000e+00</span></p></td><td class="cl-b6acb965"><p class="cl-b6acb960"><span class="cl-b6ab320c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-b6ae3dbc"><p class="cl-b6acb960"><span class="cl-b6ab320c">ti(wetness,et)</span></p></td><td class="cl-b6ae3dbd"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6acb968"><p class="cl-b6acb960"><span class="cl-b6ab320c"></span></p></td><td class="cl-b6ae3dbe"><p class="cl-b6acb961"><span class="cl-b6ab320c">436.57</span></p></td><td class="cl-b6acb96a"><p class="cl-b6acb961"><span class="cl-b6ab320c">0.000e+00</span></p></td><td class="cl-b6acb969"><p class="cl-b6acb960"><span class="cl-b6ab320c">***</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b6ae3dc4"><p class="cl-b6acb961"><span class="cl-b6ab320d">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b6ae3e93"><p class="cl-b6acb960"><span class="cl-b6ab320c"> </span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b6ae3e93"><p class="cl-b6acb960"><span class="cl-b6ab320c">(Dispersion parameter for Scaled t(3.013,0.081) family taken to be 1)</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b6ae3e93"><p class="cl-b6acb960"><span class="cl-b6ab320c">Null deviance: NULL on NULL degrees of freedom</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-b6ae3e93"><p class="cl-b6acb960"><span class="cl-b6ab320c">Residual deviance: NULL on NULL degrees of freedom</span></p></td></tr></tfoot></table></div>
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
<div class="tabwid"><style>.cl-3fe2f730{border-collapse:collapse;}.cl-3fd840d8{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-3fd840d9{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:italic;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-3fd840da{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-3fd840db{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-3fd840dc{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd840dd{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd840de{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd840df{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd840e0{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd840e1{width:61pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd840e2{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c62e{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c62f{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c630{width:61pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c631{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c632{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c633{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c634{width:61pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c635{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c636{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c637{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c638{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c84a{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c84b{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c84c{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c84d{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c84e{width:61pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c84f{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c850{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c851{width:68pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c852{width:61pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c853{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fd9c854{width:162pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-3fdb5246{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-3fe2f730'>
```
<caption class="Table Caption">(\#tab:gamresults16882)Summary of GAMcoefficients at SWQM 16882</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-3fd9c854"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd9c850"><p class="cl-3fd840da"><span class="cl-3fd840d8">Estimate</span></p></td><td class="cl-3fd9c853"><p class="cl-3fd840da"><span class="cl-3fd840d8">Standard Error</span></p></td><td class="cl-3fd9c852"><p class="cl-3fd840db"><span class="cl-3fd840d8">z value</span></p></td><td class="cl-3fd9c851"><p class="cl-3fd840db"><span class="cl-3fd840d8">Pr(&gt;|z|)</span></p></td><td class="cl-3fdb5246"><p class="cl-3fd840da"><span class="cl-3fd840d8">Signif.</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-3fd840df"><p class="cl-3fd840da"><span class="cl-3fd840d8">s(ewood_precip)</span></p></td><td class="cl-3fd840de"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840dc"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840e1"><p class="cl-3fd840db"><span class="cl-3fd840d8">7.67e+00</span></p></td><td class="cl-3fd840e0"><p class="cl-3fd840db"><span class="cl-3fd840d8">2.431e-03</span></p></td><td class="cl-3fd840dd"><p class="cl-3fd840da"><span class="cl-3fd840d8">**</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-3fd840df"><p class="cl-3fd840da"><span class="cl-3fd840d8">s(ewood_tmax)</span></p></td><td class="cl-3fd840de"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840dc"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840e1"><p class="cl-3fd840db"><span class="cl-3fd840d8">2.18e+00</span></p></td><td class="cl-3fd840e0"><p class="cl-3fd840db"><span class="cl-3fd840d8">6.979e-02</span></p></td><td class="cl-3fd840dd"><p class="cl-3fd840da"><span class="cl-3fd840d8">.</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-3fd840df"><p class="cl-3fd840da"><span class="cl-3fd840d8">s(lagPrecip)</span></p></td><td class="cl-3fd840de"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840dc"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840e1"><p class="cl-3fd840db"><span class="cl-3fd840d8">4.37e+00</span></p></td><td class="cl-3fd840e0"><p class="cl-3fd840db"><span class="cl-3fd840d8">2.071e-02</span></p></td><td class="cl-3fd840dd"><p class="cl-3fd840da"><span class="cl-3fd840d8">*</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-3fd840df"><p class="cl-3fd840da"><span class="cl-3fd840d8">s(doy)</span></p></td><td class="cl-3fd840de"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840dc"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840e1"><p class="cl-3fd840db"><span class="cl-3fd840d8">3.31e+01</span></p></td><td class="cl-3fd840e0"><p class="cl-3fd840db"><span class="cl-3fd840d8">-2.211e-06</span></p></td><td class="cl-3fd840dd"><p class="cl-3fd840da"><span class="cl-3fd840d8">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-3fd840df"><p class="cl-3fd840da"><span class="cl-3fd840d8">s(wetness)</span></p></td><td class="cl-3fd840de"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840dc"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840e1"><p class="cl-3fd840db"><span class="cl-3fd840d8">1.46e+01</span></p></td><td class="cl-3fd840e0"><p class="cl-3fd840db"><span class="cl-3fd840d8">2.385e-04</span></p></td><td class="cl-3fd840dd"><p class="cl-3fd840da"><span class="cl-3fd840d8">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-3fd840df"><p class="cl-3fd840da"><span class="cl-3fd840d8">s(et)</span></p></td><td class="cl-3fd840de"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840dc"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840e1"><p class="cl-3fd840db"><span class="cl-3fd840d8">4.50e-01</span></p></td><td class="cl-3fd840e0"><p class="cl-3fd840db"><span class="cl-3fd840d8">2.212e-01</span></p></td><td class="cl-3fd840dd"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-3fd840df"><p class="cl-3fd840da"><span class="cl-3fd840d8">ti(ewood_precip,ewood_tmax)</span></p></td><td class="cl-3fd840de"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840dc"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840e1"><p class="cl-3fd840db"><span class="cl-3fd840d8">1.89e-05</span></p></td><td class="cl-3fd840e0"><p class="cl-3fd840db"><span class="cl-3fd840d8">4.391e-01</span></p></td><td class="cl-3fd840dd"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-3fd840df"><p class="cl-3fd840da"><span class="cl-3fd840d8">ti(ewood_precip,wetness)</span></p></td><td class="cl-3fd840de"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840dc"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840e1"><p class="cl-3fd840db"><span class="cl-3fd840d8">2.20e-06</span></p></td><td class="cl-3fd840e0"><p class="cl-3fd840db"><span class="cl-3fd840d8">8.499e-01</span></p></td><td class="cl-3fd840dd"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-3fd840df"><p class="cl-3fd840da"><span class="cl-3fd840d8">ti(ewood_precip,et)</span></p></td><td class="cl-3fd840de"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840dc"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840e1"><p class="cl-3fd840db"><span class="cl-3fd840d8">5.52e+00</span></p></td><td class="cl-3fd840e0"><p class="cl-3fd840db"><span class="cl-3fd840d8">1.581e-02</span></p></td><td class="cl-3fd840dd"><p class="cl-3fd840da"><span class="cl-3fd840d8">*</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-3fd840df"><p class="cl-3fd840da"><span class="cl-3fd840d8">ti(ewood_tmax,wetness)</span></p></td><td class="cl-3fd840de"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840dc"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840e1"><p class="cl-3fd840db"><span class="cl-3fd840d8">6.19e-07</span></p></td><td class="cl-3fd840e0"><p class="cl-3fd840db"><span class="cl-3fd840d8">9.684e-01</span></p></td><td class="cl-3fd840dd"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-3fd840df"><p class="cl-3fd840da"><span class="cl-3fd840d8">ti(ewood_tmax,et)</span></p></td><td class="cl-3fd840de"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840dc"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd840e1"><p class="cl-3fd840db"><span class="cl-3fd840d8">8.37e-06</span></p></td><td class="cl-3fd840e0"><p class="cl-3fd840db"><span class="cl-3fd840d8">4.014e-01</span></p></td><td class="cl-3fd840dd"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-3fd9c632"><p class="cl-3fd840da"><span class="cl-3fd840d8">ti(wetness,et)</span></p></td><td class="cl-3fd9c631"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd9c62e"><p class="cl-3fd840da"><span class="cl-3fd840d8"></span></p></td><td class="cl-3fd9c630"><p class="cl-3fd840db"><span class="cl-3fd840d8">4.32e+00</span></p></td><td class="cl-3fd9c62f"><p class="cl-3fd840db"><span class="cl-3fd840d8">5.795e-02</span></p></td><td class="cl-3fd840e2"><p class="cl-3fd840da"><span class="cl-3fd840d8">.</span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-3fd9c638"><p class="cl-3fd840db"><span class="cl-3fd840d9">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-3fd9c84c"><p class="cl-3fd840da"><span class="cl-3fd840d8"> </span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-3fd9c84c"><p class="cl-3fd840da"><span class="cl-3fd840d8">(Dispersion parameter for binomial family taken to be 1)</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-3fd9c84c"><p class="cl-3fd840da"><span class="cl-3fd840d8">Null deviance: NULL on NULL degrees of freedom</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-3fd9c84c"><p class="cl-3fd840da"><span class="cl-3fd840d8">Residual deviance: NULL on NULL degrees of freedom</span></p></td></tr></tfoot></table></div>
```

```r
flextable::as_flextable(m4.gam) %>%
  flextable::set_caption("Summary of GAMcoefficients at SWQM 16882")
```

```{=html}
<div class="tabwid"><style>.cl-4010bb5c{border-collapse:collapse;}.cl-4006142c{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:normal;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-4006142d{font-family:'Arial';font-size:11pt;font-weight:normal;font-style:italic;text-decoration:none;color:rgba(0, 0, 0, 1.00);background-color:transparent;}.cl-4006142e{margin:0;text-align:left;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-4006142f{margin:0;text-align:right;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);padding-bottom:3pt;padding-top:3pt;padding-left:3pt;padding-right:3pt;line-height: 1;background-color:transparent;}.cl-40061430{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40061431{width:51pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40061432{width:55pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40061433{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40061434{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40061435{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40061436{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40079266{width:51pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40079267{width:55pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40079268{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40079269{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4007926a{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 0 solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4007926b{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4007926c{width:51pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4007926d{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4007926e{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4007926f{width:55pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40079270{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40079518{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40079519{width:51pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4007951a{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4007951b{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4007951c{width:55pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4007951d{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 0 solid rgba(255, 255, 255, 0.00);border-top: 0 solid rgba(255, 255, 255, 0.00);border-left: 0 solid rgba(255, 255, 255, 0.00);border-right: 0 solid rgba(255, 255, 255, 0.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4007951e{width:46pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-4007951f{width:58pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40079520{width:87pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40079521{width:144pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40079522{width:55pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}.cl-40091cf8{width:51pt;background-color:transparent;vertical-align: middle;border-bottom: 2pt solid rgba(0, 0, 0, 1.00);border-top: 2pt solid rgba(0, 0, 0, 1.00);border-left: 0 solid rgba(0, 0, 0, 1.00);border-right: 0 solid rgba(0, 0, 0, 1.00);margin-bottom:0;margin-top:0;margin-left:0;margin-right:0;}</style><table class='cl-4010bb5c'>
```
<caption class="Table Caption">(\#tab:gamresults16882)Summary of GAMcoefficients at SWQM 16882</caption>
```{=html}
<thead><tr style="overflow-wrap:break-word;"><td class="cl-40079521"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-4007951f"><p class="cl-4006142e"><span class="cl-4006142c">Estimate</span></p></td><td class="cl-40079520"><p class="cl-4006142e"><span class="cl-4006142c">Standard Error</span></p></td><td class="cl-40079522"><p class="cl-4006142f"><span class="cl-4006142c">z value</span></p></td><td class="cl-40091cf8"><p class="cl-4006142f"><span class="cl-4006142c">Pr(&gt;|z|)</span></p></td><td class="cl-4007951e"><p class="cl-4006142e"><span class="cl-4006142c">Signif.</span></p></td></tr></thead><tbody><tr style="overflow-wrap:break-word;"><td class="cl-40061435"><p class="cl-4006142e"><span class="cl-4006142c">s(doy)</span></p></td><td class="cl-40061434"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40061433"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40061432"><p class="cl-4006142f"><span class="cl-4006142c">815.535</span></p></td><td class="cl-40061431"><p class="cl-4006142f"><span class="cl-4006142c">0.0000</span></p></td><td class="cl-40061430"><p class="cl-4006142e"><span class="cl-4006142c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-40061435"><p class="cl-4006142e"><span class="cl-4006142c">s(lagPrecip)</span></p></td><td class="cl-40061434"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40061433"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40061432"><p class="cl-4006142f"><span class="cl-4006142c">65.089</span></p></td><td class="cl-40061431"><p class="cl-4006142f"><span class="cl-4006142c">0.0000</span></p></td><td class="cl-40061430"><p class="cl-4006142e"><span class="cl-4006142c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-40061435"><p class="cl-4006142e"><span class="cl-4006142c">te(ewood_precip,wetness)</span></p></td><td class="cl-40061434"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40061433"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40061432"><p class="cl-4006142f"><span class="cl-4006142c">65.459</span></p></td><td class="cl-40061431"><p class="cl-4006142f"><span class="cl-4006142c">0.0000</span></p></td><td class="cl-40061430"><p class="cl-4006142e"><span class="cl-4006142c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-40061435"><p class="cl-4006142e"><span class="cl-4006142c">te(ewood_precip,et)</span></p></td><td class="cl-40061434"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40061433"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40061432"><p class="cl-4006142f"><span class="cl-4006142c">454.991</span></p></td><td class="cl-40061431"><p class="cl-4006142f"><span class="cl-4006142c">0.0000</span></p></td><td class="cl-40061430"><p class="cl-4006142e"><span class="cl-4006142c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-40061435"><p class="cl-4006142e"><span class="cl-4006142c">te(ewood_tmax,wetness)</span></p></td><td class="cl-40061434"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40061433"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40061432"><p class="cl-4006142f"><span class="cl-4006142c">91.155</span></p></td><td class="cl-40061431"><p class="cl-4006142f"><span class="cl-4006142c">0.0000</span></p></td><td class="cl-40061430"><p class="cl-4006142e"><span class="cl-4006142c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-40061435"><p class="cl-4006142e"><span class="cl-4006142c">te(ewood_tmax,et)</span></p></td><td class="cl-40061434"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40061433"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40061432"><p class="cl-4006142f"><span class="cl-4006142c">272.787</span></p></td><td class="cl-40061431"><p class="cl-4006142f"><span class="cl-4006142c">0.0000</span></p></td><td class="cl-40061430"><p class="cl-4006142e"><span class="cl-4006142c">***</span></p></td></tr><tr style="overflow-wrap:break-word;"><td class="cl-40079268"><p class="cl-4006142e"><span class="cl-4006142c">te(wetness,et)</span></p></td><td class="cl-4007926a"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40079269"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td><td class="cl-40079267"><p class="cl-4006142f"><span class="cl-4006142c">0.251</span></p></td><td class="cl-40079266"><p class="cl-4006142f"><span class="cl-4006142c">0.1536</span></p></td><td class="cl-40061436"><p class="cl-4006142e"><span class="cl-4006142c"></span></p></td></tr></tbody><tfoot><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-40079270"><p class="cl-4006142f"><span class="cl-4006142d">Signif. codes: 0 &lt;= '***' &lt; 0.001 &lt; '**' &lt; 0.01 &lt; '*' &lt; 0.05 &lt; '.' &lt; 0.1 &lt; '' &lt; 1</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4007951a"><p class="cl-4006142e"><span class="cl-4006142c"> </span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4007951a"><p class="cl-4006142e"><span class="cl-4006142c">(Dispersion parameter for Scaled t(3,0.952) family taken to be 1)</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4007951a"><p class="cl-4006142e"><span class="cl-4006142c">Null deviance: NULL on NULL degrees of freedom</span></p></td></tr><tr style="overflow-wrap:break-word;"><td  colspan="6"class="cl-4007951a"><p class="cl-4006142e"><span class="cl-4006142c">Residual deviance: NULL on NULL degrees of freedom</span></p></td></tr></tfoot></table></div>
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
## 5 GAM                     0.601         0.373       -0.109
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
## 5 GAM                     0.636         0.452       0.402
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
## 5 GAM                      -15.8         -5.8        -12
```

## Model validation


```r
## use k-fold cross validation to estimate
## how well the approach fits to ujnknown data
## randomly split the data 80-20 to fit and validate the GAM.
library(modelr)
```

```
## 
## Attaching package: 'modelr'
```

```
## The following object is masked from 'package:gratia':
## 
##     add_residuals
```

```
## The following objects are masked from 'package:hydroGOF':
## 
##     mae, mse, rmse
```

```r
set.seed(123)

df_16396 %>%
  ## create a list of data resamples.
  ## randomly splits the dataset into 80-20
  ## n number of times
  crossv_kfold() %>%
  mutate(model = map(train, ~gam(adjusted_flow~
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
                                  data = as.data.frame(.),
                                  family = scat(link = "log"),
                                  method = "REML",
                                  control = gam.control(nthreads = 2))),
         preds = map2(test,model, ~predict(.y, newdata = as.data.frame(.x),
                                           type = "response")),
         NSE = map2_dbl(test, preds, ~hydroGOF::NSE(as.numeric(as.data.frame(.x)$adjusted_flow),
                                                    as.numeric(.y))),
         KGE = map2_dbl(test, preds, ~hydroGOF::KGE(as.numeric(as.data.frame(.x)$adjusted_flow),
                                                    as.numeric(.y))),
         MAE = map2_dbl(test, preds, ~hydroGOF::mae(as.numeric(as.data.frame(.x)$adjusted_flow),
                                                    as.numeric(.y))),
         pbias = map2_dbl(test, preds, ~hydroGOF::pbias(as.numeric(as.data.frame(.x)$adjusted_flow),
                                                    as.numeric(.y)))) -> df_16396_kfold
```

```
## Warning: Problem with `mutate()` input `model`.
## i Fitting terminated with step failure - check results carefully
## i Input `model` is `map(...)`.
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## Warning: Problem with `mutate()` input `model`.
## i Fitting terminated with step failure - check results carefully
## i Input `model` is `map(...)`.
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```
## Warning: Problem with `mutate()` input `model`.
## i Fitting terminated with step failure - check results carefully
## i Input `model` is `map(...)`.
```

```
## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
## G$L, : Fitting terminated with step failure - check results carefully
```

```r
df_16396_kfold
```

```
## # A tibble: 5 x 9
##   train        test         .id   model    preds         NSE     KGE   MAE pbias
##   <named list> <named list> <chr> <named > <named >    <dbl>   <dbl> <dbl> <dbl>
## 1 <resample [~ <resample [~ 1     <gam>    <dbl [6~ -6.84e-3 -0.489  101.  -92.6
## 2 <resample [~ <resample [~ 2     <gam>    <dbl [6~ -1.57e-2 -0.647  829.  -97.4
## 3 <resample [~ <resample [~ 3     <gam>    <dbl [6~ -2.60e-2 -0.515   22.0 -63.3
## 4 <resample [~ <resample [~ 4     <gam>    <dbl [6~ -2.13e+1 -3.27    21.0 176. 
## 5 <resample [~ <resample [~ 5     <gam>    <dbl [6~ -9.27e-1  0.0170  10.4  55.8
```

```r
# 
# m1.gam <- gam(adjusted_flow ~
#             s(ewood_precip, bs = "ts", k = 5) +
#             s(ewood_tmax, bs = "ts", k = 5) +
#             s(lagPrecip, bs = "ts", k = 5) +
#             s(doy, bs = "cc", k = 8) +
#             s(wetness, bs = "ts", k = 5) +
#             s(et, bs = "ts", k = 5) +
#               ti(ewood_precip, lagPrecip, k = 5) +
#               ti(ewood_precip, ewood_tmax, k = 5) +
#               ti(ewood_precip, wetness, k = 5) +
#               ti(ewood_precip, et, k = 5) +
#               ti(ewood_tmax, wetness, k = 5) +
#               ti(ewood_tmax, et, k = 5) +
#               ti(wetness, et, k = 5),
#           data = df_16396,
#           select = TRUE,
#           family = scat(link = "log"),
#           method = "REML",
#           control = gam.control(nthreads = 2))
```




## References
