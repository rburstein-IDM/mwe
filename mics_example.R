rm(list=ls())

data_dir <- "C:/Users/qdong/Dropbox/IDM-UW/Data/"

library(foreign)
library(data.table)
library(rgdal)
library(readr)
library(readxl)
library(raster)
library(Hmisc)
library(gdata)
library(sp)

##################################
# mics 2016 children
##################################

mics2016_chil_raw <- read.spss(paste0(data_dir, "mics_nigeria_07022018/Nigeria MICS 201617 SPSS Datasets/ch.sav"),
                               use.value.labels = T, to.data.frame = T)
# View(mics2016_chil_raw[1:10, ])

col_vt <- c("UF1", "UF2", "LN", "UF8D", "UF8M", "UF8Y", 
            "AG1D", "AG1M", "AG1Y", "AG2", 
            "IM0CA", "IM0CB", "IM0CC", "IM0CD", 
            "IM1", "IM3BD", "IM3BM", "IM3BY", "IM3IPVD", "IM3IPVM", "IM3IPVM", 
            "IM3PCV1D", "IM3PCV1M", "IM3PCV1Y", "IM3PCV2D", "IM3PCV2M", "IM3PCV2Y", "IM3PCV3D", "IM3PCV3M", "IM3PCV3Y",
            "IM30PV0D", "IM30PV0M", "IM30PV0Y", "IM30PV1D", "IM30PV1M", "IM30PV1Y", "IM30PV2D", "IM30PV2M", "IM30PV2Y", "IM30PV3D", "IM30PV3M", "IM30PV3Y", 
            "IM3PENTA1D", "IM3PENTA1M", "IM3PENTA1Y", "IM3PENTA2D", "IM3PENTA2M", "IM3PENTA2Y", "IM3PENTA3D", "IM3PENTA3M", "IM3PENTA3Y", 
            "IM3HEPB0D", "IM3HEPB0M", "IM3HEPB0Y", 
            "IM3MD", "IM3MM", "IM3MY", "IM3YD", "IM3YM", "IM3YY", 
            "IM3V1D", "IM3V1M", "IM3V1Y", "IM3V2D", "IM3V2M", "IM3V2Y", 
            "IM8", "IM10", "IM10A", "IM12A", "IM12B", "IM15A", "IM15B", "IM16", "IM17",
            "HH6", "HH7", "chweight", "Zone", "CDOI", "CDOB", "CAGE", "CAGED", "chweightkano", "chweightlagos")

col_name_vt <- c("cluster", "hh", "child", "interv_d", "interv_m", "interv_y",
                 "birth_d", "birth_m", "birth_y", "age_y",
                 "camp_a", "camp_b", "camp_c", "camp_d",
                 "hascard", "bcg_d", "bcg_m", "bcg_y", "ipv_d", "ipv_m", "ipv_y",
                 "pcv1_d", "pcv1_m", "pcv1_y", "pcv2_d", "pcv2_m", "pcv2_y", "pcv3_d", "pcv3_m", "pcv3_y", 
                 "polio0_d", "polio0_m", "polio0_y", "polio1_d", "polio1_m", "polio1_y", "polio2_d", "polio2_m", "polio2_y", "polio3_d", "polio3_m", "polio3_y",
                 "dpt1_d", "dpt1_m", "dpt1_y", "dpt2_d", "dpt2_m", "dpt2_y", "dpt3_d", "dpt3_m", "dpt3_y", 
                 "hb0_d", "hb0_m", "hb0_y", 
                 "measles_d", "measles_m", "measles_y", "yf_d", "yf_m", "yf_y", 
                 "va1_d", "va1_m", "va1_y", "va2_d", "va2_m", "va2_y", 
                 "polio", "polio_num", "ipv", "dpt", "dpt_num", "pcv", "pcv_num", "measles", "yf",
                 "type", "state", "chsampwgt", "region", "interv_date", "birth_date", "age_m", "age_d", "chsampwgt_kano", "chsampwgt_lagos")

mics2016_chil_dt <- as.data.table(mics2016_chil_raw[, col_vt])
setnames(mics2016_chil_dt, col_name_vt)

mics2016_chil_dt[, "DPT3" := ifelse(!is.na(dpt) & dpt == "Yes" & !dpt_num %in% c("1", "2", "Missing"), 1, 
                                    ifelse(!is.na(dpt3_d) & !dpt3_d %in% c("Not given", "DK"), 1, 0))]

