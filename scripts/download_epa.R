library(dplyr)
library(echor)
library(ggplot2)
library(tsibble)
library(units)

## swqm 16396 includes TX0025071 (Still Creek WWTF) and TX0113603 (Sanderson)
## swqm 16882 includes TX0025071 (Still Creek WWTF)

stillcreek <- echoGetEffluent(p_id = "TX0025071",
                              parameter_code = '50050',
                              start_date = "02/01/2020",
                              end_date = "03/01/2021")

sanderson <- echoGetEffluent(p_id = "TX0113603",
                             parameter_code = '50050',
                             start_date = "02/01/2020",
                             end_date = "03/01/2021")

stillcreek %>%
  bind_rows(sanderson) %>%
  mutate(monitoring_period_end_date = as.Date(monitoring_period_end_date,
                                              format = "%m/%d/%Y"),
         dmr_value_nmbr = as.numeric(dmr_value_nmbr)) %>%
  filter(statistical_base_code == "DB") %>%
  ggplot() +
  geom_step(aes(monitoring_period_end_date, dmr_value_nmbr, color = npdes_id))


stillcreek %>%
  bind_rows(sanderson) %>%
  mutate(monitoring_period_end_date = as.Date(monitoring_period_end_date,
                                              format = "%m/%d/%Y"),
         dmr_value_nmbr = as.numeric(dmr_value_nmbr)) %>%
  filter(statistical_base_code == "DB") %>%
  select(npdes_id, perm_feature_nmbr, monitoring_period_end_date, dmr_value_nmbr) %>%
  as_tsibble(key = npdes_id, index = monitoring_period_end_date) %>%
  group_by_key() %>%
  fill_gaps() %>%
  tidyr::fill(dmr_value_nmbr, .direction = "up") %>%
  as_tibble() %>%
  select(npdes_id, date = monitoring_period_end_date, mgd = dmr_value_nmbr) %>%
  mutate(mgd = set_units(mgd, "1E6 * gallon/day")) %>%
  mutate(cfs = set_units(mgd, "foot^3/s")) %>%
  readr::write_csv("Data/EPA_WWTF/mean_daily_discharges.csv")
