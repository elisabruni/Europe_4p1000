module load cdo
#Calculate average npp across models
cdo ensmean /RCP26/DAILY_FORCING/pr_day_2006_2100.nc4 /RCP26/DAILY_FORCING/MIROC/pr_day_MIROC_rcp26_2006_2100.nc4 pr_day_ORCH_MIROC_rcp26_2006_2100.nc4

cdo ensmean /RCP26/DAILY_FORCING/prsn_day_2006_2100.nc4 /RCP26/DAILY_FORCING/MIROC/prsn_day_MIROC_rcp26_2006_2100.nc4 prsn_day_ORCH_MIROC_rcp26_2006_2100.nc4

cdo ensmean /RCP26/DAILY_FORCING/tas_day_2006_2100.nc4 /RCP26/DAILY_FORCING/MIROC/tas_day_MIROC_rcp26_2006_2100.nc4 tas_day_ORCH_MIROC_rcp26_2006_2100.nc4
