#################
# Preparation of data for DHIN analysis
# Citation: Goldstein ND, Zhang Y, Burstyn I, Vasireddy K, Garland D, Enright M, Siegel SD. Towards improved methodology to measure the uptake of mammography in populations using insurance claims: a case study from Delaware, USA. Manuscript in preparation.
# 9/30/25 -- Neal Goldstein (adapted from code from Yuchen Zhang)
#################


#functions
library(sf)
library(dplyr)
library(stringr)
library(tidycensus)
library(tidyverse)


path = "<path to files>"
path_datasource = "<path to files>"

#census tract boundaries
DE_2020_tracts<-st_read(paste0(path_datasource,"DE_tract_2020/tl_2020_10_tract.shp"))
DE_2010_tracts<-st_read(paste0(path_datasource,"DE_DEStatePlaneProjection/DE_Tract_2010_SP.shp"))

#mammo data
DHIN_mammo_data = read.csv(paste0(path_datasource,"/DHIN_data_aggregated.csv"), as.is=T)
DHIN_mammo_data$GEOID = as.character(DHIN_mammo_data$GEOID)

##alghough the dataframe name below are NCC, they actually reprensents the whole DE state.
NCC_tract_decennial_2010<-get_decennial(#--------
                                        geography = "tract",  #at census tract level,  this can also be state, or other level
                                        variables =  c(
                                          total_female_35_39='P012037',
                                          total_female_40_44='P012038',
                                          total_female_45_49='P012039',
                                          total_female_50_54='P012040',
                                          total_female_55_59='P012041',
                                          total_female_60_61='P012042',
                                          total_female_62_64='P012043',
                                          total_female_65_66='P012044',
                                          total_female_67_69='P012045',
                                          total_female_70_74='P012046',
                                          total_female_75_79='P012047',
                                          total_female_80_84='P012048',
                                          total_female_85andover='P012049'
                                          
                                        ), 
                                        state= 10,   #Delaware FIPS
                                        county=c(001,003,005),  # NCC FIPS
                                        year = 2010, 
                                        geometry=TRUE
)

#aggregate block group level census data:
NCC_tract_2010_aggregated<- NCC_tract_decennial_2010%>%spread(variable, value) %>%filter(GEOID != '100039901000') # remove the tract with zero population

#clean the dataset, create aggregated values:
NCC_tract_2010_aggregated<-NCC_tract_2010_aggregated%>%
  mutate(
    total_female_over40=total_female_40_44+total_female_45_49+
      total_female_50_54+total_female_55_59+total_female_60_61+
      total_female_62_64+total_female_65_66+total_female_67_69+
      total_female_70_74+total_female_75_79+total_female_80_84+
      total_female_85andover
  )%>%
  select(
    GEOID, total_female_over40
  )


#aggregate tract level census data:-----------------
NCC_tract_2010_aggregated<- NCC_tract_decennial_2010%>%spread(variable, value) %>%filter(GEOID != '100039901000') # remove the tract with zero population

#clean the dataset, create aggregated values:
NCC_tract_2010_aggregated<-NCC_tract_2010_aggregated%>%
  mutate(
    total_female_over40=total_female_40_44+total_female_45_49+
      total_female_50_54+total_female_55_59+total_female_60_61+
      total_female_62_64+total_female_65_66+total_female_67_69+
      total_female_70_74+total_female_75_79+total_female_80_84+
      total_female_85andover,
    total_female_40to44_portion_of_35to44=round (total_female_40_44/(total_female_35_39+total_female_40_44),digits=2)
  )%>%
  select(
    GEOID, total_female_over40, total_female_40to44_portion_of_35to44
  )

total_women_age40andup_census<-sum(NCC_tract_2010_aggregated$total_female_over40)


# bring in 2020 census data: ----------------------------------------------
NCC_tract_decennial_2020<-get_decennial(#------
                                        geography = "tract",  #at census tract level,  this can also be state, or other level
                                        variables =  c(
                                          
                                          # add age 35+39 to be able to caluclate the proportion of 40-44 patients for public insurance coverage data.
                                          total_female_35_39_white='P12A_037N',
                                          total_female_35_39_black='P12B_037N',
                                          total_female_35_39_americanindianalskan='P12C_037N',
                                          total_female_35_39_asian='P12D_037N',
                                          total_female_35_39_nativehawaii='P12E_037N',
                                          total_female_35_39_otherrace='P12F_037N',
                                          total_female_35_39_twomorerace='P12G_037N',
                                          
                                          
                                          #--- female by age and race:  
                                          total_female_40_44_white='P12A_038N',
                                          total_female_40_44_black='P12B_038N',
                                          total_female_40_44_americanindianalskan='P12C_038N',
                                          total_female_40_44_asian='P12D_038N',
                                          total_female_40_44_nativehawaii='P12E_038N',
                                          total_female_40_44_otherrace='P12F_038N',
                                          total_female_40_44_twomorerace='P12G_038N',
                                          
                                          total_female_45_49_white='P12A_039N',
                                          total_female_45_49_black='P12B_039N',
                                          total_female_45_49_americanindianalskan='P12C_039N',
                                          total_female_45_49_asian='P12D_039N',
                                          total_female_45_49_nativehawaii='P12E_039N',
                                          total_female_45_49_otherrace='P12F_039N',
                                          total_female_45_49_twomorerace='P12G_039N',
                                          
                                          total_female_50_54_white='P12A_040N',
                                          total_female_50_54_black='P12B_040N',
                                          total_female_50_54_americanindianalskan='P12C_040N',
                                          total_female_50_54_asian='P12D_040N',
                                          total_female_50_54_nativehawaii='P12E_040N',
                                          total_female_50_54_otherrace='P12F_040N',
                                          total_female_50_54_twomorerace='P12G_040N',
                                          
                                          total_female_55_59_white='P12A_041N',
                                          total_female_55_59_black='P12B_041N',
                                          total_female_55_59_americanindianalskan='P12C_041N',
                                          total_female_55_59_asian='P12D_041N',
                                          total_female_55_59_nativehawaii='P12E_041N',
                                          total_female_55_59_otherrace='P12F_041N',
                                          total_female_55_59_twomorerace='P12G_041N',
                                          
                                          total_female_60_61_white='P12A_042N',
                                          total_female_60_61_black='P12B_042N',
                                          total_female_60_61_americanindianalskan='P12C_042N',
                                          total_female_60_61_asian='P12D_042N',
                                          total_female_60_61_nativehawaii='P12E_042N',
                                          total_female_60_61_otherrace='P12F_042N',
                                          total_female_60_61_twomorerace='P12G_042N',
                                          
                                          total_female_62_64_white='P12A_043N',
                                          total_female_62_64_black='P12B_043N',
                                          total_female_62_64_americanindianalskan='P12C_043N',
                                          total_female_62_64_asian='P12D_043N',
                                          total_female_62_64_nativehawaii='P12E_043N',
                                          total_female_62_64_otherrace='P12F_043N',
                                          total_female_62_64_twomorerace='P12G_043N',
                                          
                                          total_female_65_66_white='P12A_044N',
                                          total_female_65_66_black='P12B_044N',
                                          total_female_65_66_americanindianalskan='P12C_044N',
                                          total_female_65_66_asian='P12D_044N',
                                          total_female_65_66_nativehawaii='P12E_044N',
                                          total_female_65_66_otherrace='P12F_044N',
                                          total_female_65_66_twomorerace='P12G_044N',
                                          
                                          total_female_67_69_white='P12A_045N',
                                          total_female_67_69_black='P12B_045N',
                                          total_female_67_69_americanindianalskan='P12C_045N',
                                          total_female_67_69_asian='P12D_045N',
                                          total_female_67_69_nativehawaii='P12E_045N',
                                          total_female_67_69_otherrace='P12F_045N',
                                          total_female_67_69_twomorerace='P12G_045N',
                                          
                                          total_female_70_74_white='P12A_046N',
                                          total_female_70_74_black='P12B_046N',
                                          total_female_70_74_americanindianalskan='P12C_046N',
                                          total_female_70_74_asian='P12D_046N',
                                          total_female_70_74_nativehawaii='P12E_046N',
                                          total_female_70_74_otherrace='P12F_046N',
                                          total_female_70_74_twomorerace='P12G_046N',
                                          
                                          total_female_75_79_white='P12A_047N',
                                          total_female_75_79_black='P12B_047N',
                                          total_female_75_79_americanindianalskan='P12C_047N',
                                          total_female_75_79_asian='P12D_047N',
                                          total_female_75_79_nativehawaii='P12E_047N',
                                          total_female_75_79_otherrace='P12F_047N',
                                          total_female_75_79_twomorerace='P12G_047N',
                                          
                                          total_female_80_84_white='P12A_048N',
                                          total_female_80_84_black='P12B_048N',
                                          total_female_80_84_americanindianalskan='P12C_048N',
                                          total_female_80_84_asian='P12D_048N',
                                          total_female_80_84_nativehawaii='P12E_048N',
                                          total_female_80_84_otherrace='P12F_048N',
                                          total_female_80_84_twomorerace='P12G_048N',
                                          
                                          total_female_85andover_white='P12A_049N',
                                          total_female_85andover_black='P12B_049N',
                                          total_female_85andover_americanindianalskan='P12C_049N',
                                          total_female_85andover_asian='P12D_049N',
                                          total_female_85andover_nativehawaii='P12E_049N',
                                          total_female_85andover_otherrace='P12F_049N',
                                          total_female_85andover_twomorerace='P12G_049N'
                                          #--- female by age and race end 
                                        ), 
                                        state= 10,   #Delaware FIPS
                                        county=c(001,003,005),  # NCC FIPS
                                        year = 2020, #only work for 2010 not 2020. is 2020 not available maybe?
                                        sumfile = "dhc",  # to validate against dhc file
                                        geometry=TRUE
)

#get to know 2020 variables: 
decennial_2020_vars <- load_variables(
  year = 2020, 
  dataset = "dhc",#demographics and housing characteristics   #"pl", 
  cache = FALSE
)


#validation of these numbers via the GUI census : https://data.census.gov/table/DECENNIALDHC2020.P12?g=050XX00US10003$1400000&y=2020&d=DEC%20Demographic%20and%20Housing%20Characteristics

#aggregate tract  level census data:--------
NCC_tract_2020_aggregated<- NCC_tract_decennial_2020%>%spread(variable, value) %>%filter(GEOID != '100039901000') # remove the tract with zero population

#clean the dataset, create aggregated values:
NCC_tract_2020_aggregated<-NCC_tract_2020_aggregated%>%
  mutate(
    total_female_35_39= total_female_35_39_white+
      total_female_35_39_black+
      total_female_35_39_americanindianalskan+
      total_female_35_39_asian+
      total_female_35_39_nativehawaii+
      total_female_35_39_otherrace+
      total_female_35_39_twomorerace,
    
    
    total_female_40_44= total_female_40_44_white+
      total_female_40_44_black+
      total_female_40_44_americanindianalskan+
      total_female_40_44_asian+
      total_female_40_44_nativehawaii+
      total_female_40_44_otherrace+
      total_female_40_44_twomorerace,
    
    total_female_45_49= total_female_45_49_white+
      total_female_45_49_black+
      total_female_45_49_americanindianalskan+
      total_female_45_49_asian+
      total_female_45_49_nativehawaii+
      total_female_45_49_otherrace+
      total_female_45_49_twomorerace,
    
    total_female_50_54= total_female_50_54_white+
      total_female_50_54_black+
      total_female_50_54_americanindianalskan+
      total_female_50_54_asian+
      total_female_50_54_nativehawaii+
      total_female_50_54_otherrace+
      total_female_50_54_twomorerace,
    
    total_female_55_59= total_female_55_59_white+
      total_female_55_59_black+
      total_female_55_59_americanindianalskan+
      total_female_55_59_asian+
      total_female_55_59_nativehawaii+
      total_female_55_59_otherrace+
      total_female_55_59_twomorerace,
    
    total_female_60_61= total_female_60_61_white+
      total_female_60_61_black+
      total_female_60_61_americanindianalskan+
      total_female_60_61_asian+
      total_female_60_61_nativehawaii+
      total_female_60_61_otherrace+
      total_female_60_61_twomorerace,
    
    total_female_62_64= total_female_62_64_white+
      total_female_62_64_black+
      total_female_62_64_americanindianalskan+
      total_female_62_64_asian+
      total_female_62_64_nativehawaii+
      total_female_62_64_otherrace+
      total_female_62_64_twomorerace,
    
    total_female_65_66= total_female_65_66_white+
      total_female_65_66_black+
      total_female_65_66_americanindianalskan+
      total_female_65_66_asian+
      total_female_65_66_nativehawaii+
      total_female_65_66_otherrace+
      total_female_65_66_twomorerace,
    
    total_female_67_69= total_female_67_69_white+
      total_female_67_69_black+
      total_female_67_69_americanindianalskan+
      total_female_67_69_asian+
      total_female_67_69_nativehawaii+
      total_female_67_69_otherrace+
      total_female_67_69_twomorerace,
    
    total_female_70_74= total_female_70_74_white+
      total_female_70_74_black+
      total_female_70_74_americanindianalskan+
      total_female_70_74_asian+
      total_female_70_74_nativehawaii+
      total_female_70_74_otherrace+
      total_female_70_74_twomorerace,
    
    total_female_75_79= total_female_75_79_white+
      total_female_75_79_black+
      total_female_75_79_americanindianalskan+
      total_female_75_79_asian+
      total_female_75_79_nativehawaii+
      total_female_75_79_otherrace+
      total_female_75_79_twomorerace,
    
    total_female_80_84= total_female_80_84_white+
      total_female_80_84_black+
      total_female_80_84_americanindianalskan+
      total_female_80_84_asian+
      total_female_80_84_nativehawaii+
      total_female_80_84_otherrace+
      total_female_80_84_twomorerace,
    
    total_female_85andover= total_female_85andover_white+
      total_female_85andover_black+
      total_female_85andover_americanindianalskan+
      total_female_85andover_asian+
      total_female_85andover_nativehawaii+
      total_female_85andover_otherrace+
      total_female_85andover_twomorerace
    
  )%>%
  mutate(
    total_female_over40=total_female_40_44+total_female_45_49+total_female_50_54+
      total_female_55_59+total_female_60_61+total_female_62_64+total_female_65_66+
      total_female_67_69+total_female_70_74+total_female_75_79+total_female_80_84+
      total_female_85andover,
    total_female_40to44_portion_of_35to44=round (total_female_40_44/(total_female_35_39+total_female_40_44),digits=2)
    
  )%>%
  # select(GEOID, NAME, total_female_over40, total_female_40_44,total_female_45_49,
  #        total_female_50_54,total_female_55_59,total_female_60_61,total_female_62_64,
  #        total_female_65_66,total_female_67_69,total_female_70_74,total_female_75_79,
  #        total_female_80_84,total_female_85andover)
  select(GEOID,total_female_over40,total_female_40to44_portion_of_35to44 )


#for 2019 mammo rate, use 2010 census districting and census------
mammo_data_2019<-NCC_tract_2010_aggregated%>%left_join( DHIN_mammo_data%>%filter(year==2019), by='GEOID'  )%>%
  mutate(
    year= ifelse(is.na(year), 2019, year)
  )

#for 2022 mammo rate, use 2020 census districting and census------
mammo_data_2022<-NCC_tract_2020_aggregated%>%left_join( DHIN_mammo_data%>%filter(year==2022), by='GEOID'  )

#explore distribution of DHIN membership
hist(mammo_data_2019$total_members/mammo_data_2019$total_female_over40, breaks="fd")
hist(mammo_data_2022$total_members/mammo_data_2022$total_female_over40, breaks="fd")
hist(c(mammo_data_2019$total_members/mammo_data_2019$total_female_over40, mammo_data_2022$total_members/mammo_data_2022$total_female_over40), breaks="fd")
median(mammo_data_2019$total_members/mammo_data_2019$total_female_over40, na.rm=T)
median(mammo_data_2022$total_members/mammo_data_2022$total_female_over40, na.rm=T)
median(c(mammo_data_2019$total_members/mammo_data_2019$total_female_over40, mammo_data_2022$total_members/mammo_data_2022$total_female_over40), na.rm=T)
hist(mammo_data_2019$total_members[mammo_data_2019$total_members<mammo_data_2019$total_female_over40]/mammo_data_2019$total_female_over40[mammo_data_2019$total_members<mammo_data_2019$total_female_over40], breaks="fd")
hist(mammo_data_2022$total_members[mammo_data_2022$total_members<mammo_data_2022$total_female_over40]/mammo_data_2022$total_female_over40[mammo_data_2022$total_members<mammo_data_2022$total_female_over40], breaks="fd")
hist(c(mammo_data_2019$total_members[mammo_data_2019$total_members<mammo_data_2019$total_female_over40]/mammo_data_2019$total_female_over40[mammo_data_2019$total_members<mammo_data_2019$total_female_over40], mammo_data_2022$total_members[mammo_data_2022$total_members<mammo_data_2022$total_female_over40]/mammo_data_2022$total_female_over40[mammo_data_2022$total_members<mammo_data_2022$total_female_over40]), breaks="fd")
median(mammo_data_2019$total_members[mammo_data_2019$total_members<mammo_data_2019$total_female_over40]/mammo_data_2019$total_female_over40[mammo_data_2019$total_members<mammo_data_2019$total_female_over40], na.rm=T)
median(mammo_data_2022$total_members[mammo_data_2022$total_members<mammo_data_2022$total_female_over40]/mammo_data_2022$total_female_over40[mammo_data_2022$total_members<mammo_data_2022$total_female_over40], na.rm=T)
median(c(mammo_data_2019$total_members[mammo_data_2019$total_members<mammo_data_2019$total_female_over40]/mammo_data_2019$total_female_over40[mammo_data_2019$total_members<mammo_data_2019$total_female_over40], mammo_data_2022$total_members[mammo_data_2022$total_members<mammo_data_2022$total_female_over40]/mammo_data_2022$total_female_over40[mammo_data_2022$total_members<mammo_data_2022$total_female_over40]), na.rm=T)
cap_rate=1 #cap at 100%, or perfect sensitivity (in other words, do not adjust capped tracts)


#create final data sets for mapping
final_data_geo<-rbind(mammo_data_2019, mammo_data_2022)%>%
  mutate(
    DHIN_total_members=if_else( total_members>= total_female_over40,as.integer(total_female_over40*cap_rate), total_members ),
    cap_status=  if_else( total_members>= total_female_over40,'capped', 'not capped' )
  )%>%
  mutate(
    
    DHIN_pop_proportion=round(DHIN_total_members/total_female_over40,digits=2),#used capped DHIN total members, instead of raw DHIN "total_members"
    DHIN_mammo_census_rate=round(mammo_members/total_female_over40, digits=2),
    DHIN_mammo_screening_rate=round(mammo_members/DHIN_total_members, digits=2),  #used capped DHIN total members, instead of raw DHIN "total_members"
    Mammo_per_1000_census=mammo_members/(total_female_over40/1000)
    
  )%>%
  select(
    GEOID, CensusFemaleOver40=total_female_over40, 
    total_female_40to44_portion_of_35to44,
    year,
    DHIN_total_members,
    DHIN_pop_proportion,
    cap_status,
    #address_estimated_members, ### NDG CODE ###
    #approximation_percentage, ### NDG CODE ###
    DHIN_mammo_members=mammo_members,
    Mammo_per_1000_census,
    DHIN_mammo_screening_rate,
    DHIN_mammo_census_rate
  ) 


final_data<-final_data_geo%>%st_drop_geometry() %>%filter(!is.na(DHIN_total_members)  )  # some 2020 census tracts didn't exist for 2019 DHIN data, thus they are NA values, remove them as they are not relavent for 2022 analysis


#add census variables needed for simulation routine

# census public coverage insurance pct: -------
#2022
# 
# ACS_5years_vars <- load_variables(
#   year = 2022, 
#   dataset = "acs5",#demographics and housing characteristics   #"pl", 
#   cache = FALSE
# )%>%filter(concept=='Public Health Insurance Status by Sex by Age' & geography=='tract')


ACS_five_year_gov_health_insurance_raw_2022<-get_acs(
  geography="tract",
  #table="B27003",
  variables= c(
    "B27003_044", #female 35-44 wt public insurance coverage
    "B27003_047", #45-54
    "B27003_050", #55-64
    "B27003_053", #65-74
    "B27003_056"  #75 and over
  ),
  cache_table = TRUE,
  year = 2022,
  #output = "tidy",
  output = "wide",
  state = 10,
  geometry = FALSE,
  survey = "acs5"
  
)

# census_gov_health_insurance_raw_2022<-get_acs(
#   geography="tract",
#   variables = c(
#     "DP03_0095", #Civilian non institutionalized population total
#     "DP03_0096", #under 0095, the ones with health insurance coverage(including private and public coverage)
#     "DP03_0099", #under 0095, the ones with NO health insurance coverage, 0096+0099 should equal to 0095
#     
#     "DP03_0097", #with private insurance
#     "DP03_0098", # with public insurance,  for some reason, 0097+0098 is not the same as 0096.  
#     
#     "DP03_0098P" # the actual variable we need, the percentage of public insurance coverage
#   ),
#   #table = "DP03",  #once variables are specified, table can't be specified anymore. 
#   cache_table = TRUE,
#   year = 2022,
#   #output = "tidy",
#   output = "wide",
#   state = 10,
#   geometry = FALSE,
#   survey = "acs5"
# )
census_gov_health_insurance_2022<-ACS_five_year_gov_health_insurance_raw_2022%>%mutate(
  female_wt_public_coverage_35_44=B27003_044E, #female 35-44 wt public insuarance coverage
  female_wt_public_coverage_45_54=B27003_047E, #female 45-54 wt public insuarance coverage
  female_wt_public_coverage_55_64=B27003_050E, #female 55-64 wt public insuarance coverage#
  female_wt_public_coverage_65_74=B27003_053E, #female 65-74vwt public insuarance coverage#
  female_wt_public_coverage_75andover=B27003_056E  #female 75 and over wt public insuarance coverage#
)%>%select(GEOID, NAME,female_wt_public_coverage_35_44,female_wt_public_coverage_45_54,
           female_wt_public_coverage_55_64,female_wt_public_coverage_65_74,
           female_wt_public_coverage_75andover)%>%  mutate_all(., ~replace(., is.na(.), 0))%>%
  left_join(NCC_tract_2020_aggregated, by="GEOID")%>%
  mutate(
    female_wt_public_coverage_40_44 =as.integer( female_wt_public_coverage_35_44*total_female_40to44_portion_of_35to44)
  )%>%mutate(
    FemaleOver40_wt_public_coverage=female_wt_public_coverage_40_44+female_wt_public_coverage_45_54+
      female_wt_public_coverage_55_64+female_wt_public_coverage_65_74+female_wt_public_coverage_75andover,
    ACSYear='2018-2022'
  )%>%dplyr::select(
    GEOID, ACSYear,FemaleOver40_wt_public_coverage#, total_female_over40
  )


#2019
# census_gov_health_insurance_raw_2019<-get_acs(
#   geography="tract",
#   variables = c(
#     "DP03_0095", #Civilian non institutionalized population total
#     "DP03_0096", #under 0095, the ones with health insurance coverage(including private and public coverage)
#     "DP03_0099", #under 0095, the ones with NO health insurance coverage, 0096+0099 should equal to 0095
#     
#     "DP03_0097", #with private insurance
#     "DP03_0098", # with public insurance,  for some reason, 0097+0098 is not the same as 0096.  
#     
#     "DP03_0098P" # the actual variable we need, the percentage of public insurance coverage
#   ),
#   #table = "DP03",  #once variables are specified, table can't be specified anymore. 
#   cache_table = TRUE,
#   year = 2019,
#   #output = "tidy",
#   output = "wide",
#   state = 10,
#   geometry = FALSE,
#   survey = "acs5"
# )
# #ad hoc, get county level percentage:
# census_gov_health_insurance_NCC_county_2022<-get_acs(
#   geography="county",
#   variables = c(
#     "DP03_0098P" # the actual vairable we need, the percentage of public insurance coverage
#   ),
#   #table = "DP03",  #once variables are specified, table can't be specified anymore. 
#   cache_table = TRUE,
#   year = 2022,
#   #output = "tidy",
#   output = "wide",
#   state = 10,
#   county=003,
#   geometry = FALSE,
#   survey = "acs5"
# )


ACS_five_year_gov_health_insurance_raw_2019<-get_acs(
  geography="tract",
  #table="B27003",
  variables= c(
    "B27003_044", #female 35-44 wt public insurance coverage
    "B27003_047", #45-54
    "B27003_050", #55-64
    "B27003_053", #65-74
    "B27003_056"  #75 and over
  ),
  cache_table = TRUE,
  year = 2019,
  #output = "tidy",
  output = "wide",
  state = 10,
  geometry = FALSE,
  survey = "acs5"
  
)
census_gov_health_insurance_2019<-ACS_five_year_gov_health_insurance_raw_2019%>%mutate(
  female_wt_public_coverage_35_44=B27003_044E, #female 35-44 wt public insuarance coverage
  female_wt_public_coverage_45_54=B27003_047E, #female 45-54 wt public insuarance coverage
  female_wt_public_coverage_55_64=B27003_050E, #female 55-64 wt public insuarance coverage#
  female_wt_public_coverage_65_74=B27003_053E, #female 65-74vwt public insuarance coverage#
  female_wt_public_coverage_75andover=B27003_056E  #female 75 and over wt public insuarance coverage#
)%>%select(GEOID, NAME,female_wt_public_coverage_35_44,female_wt_public_coverage_45_54,
           female_wt_public_coverage_55_64,female_wt_public_coverage_65_74,
           female_wt_public_coverage_75andover)%>%  mutate_all(., ~replace(., is.na(.), 0))%>%
  left_join(NCC_tract_2010_aggregated, by="GEOID")%>%
  mutate(
    female_wt_public_coverage_40_44 =as.integer( female_wt_public_coverage_35_44*total_female_40to44_portion_of_35to44)
  )%>%mutate(
    FemaleOver40_wt_public_coverage=female_wt_public_coverage_40_44+female_wt_public_coverage_45_54+
      female_wt_public_coverage_55_64+female_wt_public_coverage_65_74+female_wt_public_coverage_75andover,
    ACSYear='2015-2019'
  )%>%dplyr::select(
    GEOID, ACSYear,FemaleOver40_wt_public_coverage #, total_female_over40
  )

#merge the census data to main dataset:----------

year2022<-final_data_geo%>%filter(year=='2022')
year2022<-year2022%>%left_join(census_gov_health_insurance_2022, by="GEOID")%>%
  mutate(
    DHIN_representation_rate=round(DHIN_total_members/CensusFemaleOver40, digits = 2),
    PublicInsuranceProportion= round(FemaleOver40_wt_public_coverage/CensusFemaleOver40, digits = 2)
    
  )%>%
  mutate(
    DHIN_representation_rate_category=case_when(
      DHIN_representation_rate<0.2 ~"<0.2",
      DHIN_representation_rate<0.4 ~ "0.2-0.4",
      DHIN_representation_rate<0.6 ~ "0.4-0.6",
      DHIN_representation_rate<0.8 ~ "0.6-0.8",
      DHIN_representation_rate<1.0 ~ "0.8-1.0",
      TRUE ~ 'No data'
    ),
    
    PublicInsuranceProportion_category=case_when(
      PublicInsuranceProportion<0.2 ~"<0.2",
      PublicInsuranceProportion<0.4 ~ "0.2-0.4",
      PublicInsuranceProportion<0.6 ~ "0.4-0.6",
      PublicInsuranceProportion<0.8 ~ "0.6-0.8",
      PublicInsuranceProportion<1.0 ~ "0.8-1.0",
      TRUE ~ 'No data'
    )
  )
#map the data out, check area with higher public coverage:
# ggplot()+
#   geom_sf(data = year2022, aes( fill=PublicInsuranceProportion)) + 
#   scale_fill_viridis_c(direction=-1)+   # color shows from lightest to darkest, if 1(defualt), it's from darkest to lightest
#   xlab("Easting")+
#   ylab("Northing")+
#   ggtitle(paste0("Proportion of public insurance coverage by tract ","\n",
#                  "Reference: ACS 5 year 2022"))
# % Public insurance coverage:
# ggplot()+
#   geom_sf(data = year2022, aes( fill= PublicInsuranceProportion_category)) + 
#   # scale_fill_viridis_c(direction=-1)+   # color shows from lightest to darkest, if 1(defualt), it's from darkest to lightest
#   scale_fill_manual(values = c("<0.2" = "#7B886B", "0.2-0.4" = "#88BB92", "0.4-0.6" = "#94DDBC","0.6-0.8"="#A0ECD0", "0.8-1.0"="white",
#                                "No data"="black"))+ # Assign specific colors to categories
# 
#   xlab("Easting")+
#   ylab("Northing")+
#   ggtitle(paste0("Proportion of public insurance coverage by tract ","\n",
#                  "Reference: ACS 5 year 2022"))
# 
# 
# #% DHIN representation
# ggplot()+
#   geom_sf(data = year2022, aes( fill= DHIN_representation_rate_category)) + 
#   # scale_fill_viridis_c(direction=-1)+   # color shows from lightest to darkest, if 1(defualt), it's from darkest to lightest
#   scale_fill_manual(values = c("<0.2" = "#7B886B", "0.2-0.4" = "#88BB92", "0.4-0.6" = "#94DDBC","0.6-0.8"="#A0ECD0","0.8-1.0"="white",
#                                "No data"="black"))+ # Assign specific colors to categories
#   
#   xlab("Easting")+
#   ylab("Northing")+
#   ggtitle(paste0("Proportion of DHIN representation by tract ","\n",
#                  "Reference: ACS 5 year 2022"))


#gather dataset 2019
year2019<-final_data_geo%>%filter(year=='2019')
year2019<-year2019%>%left_join(census_gov_health_insurance_2019, by="GEOID")%>%
  mutate(
    DHIN_representation_rate=round(DHIN_total_members/CensusFemaleOver40, digits = 2),
    PublicInsuranceProportion= round(FemaleOver40_wt_public_coverage/CensusFemaleOver40, digits = 2)
    
  )%>%mutate(
    DHIN_representation_rate_category=case_when(
      DHIN_representation_rate<0.2 ~"<0.2",
      DHIN_representation_rate<0.4 ~ "0.2-0.4",
      DHIN_representation_rate<0.6 ~ "0.4-0.6",
      DHIN_representation_rate<0.8 ~ "0.6-0.8",
      DHIN_representation_rate<1.0 ~ "0.8-1.0",
      TRUE ~ 'No data'
    ),
    # #the category below has finer brackets, to show variation on its own:
    # DHIN_representation_rate_category_variation=case_when(
    #   DHIN_representation_rate<0.1 ~"<0.1",
    #   DHIN_representation_rate<0.2 ~"0.1-0.2",  
    #   DHIN_representation_rate<0.3 ~"0.2-0.3",
    #   DHIN_representation_rate<0.4 ~ "0.3-0.4",
    #   DHIN_representation_rate<0.6 ~ "0.4-0.6",
    #   DHIN_representation_rate<0.8 ~ "0.6-0.8",
    #   DHIN_representation_rate<1.0 ~ "0.8-1.0",
    #   TRUE ~ 'No data'
    # ),
    
    PublicInsuranceProportion_category=case_when(
      PublicInsuranceProportion<0.2 ~"<0.2",
      PublicInsuranceProportion<0.4 ~ "0.2-0.4",
      PublicInsuranceProportion<0.6 ~ "0.4-0.6",
      PublicInsuranceProportion<0.8 ~ "0.6-0.8",
      PublicInsuranceProportion<1.0 ~ "0.8-1.0",
      TRUE ~ 'No data'
    )
    # #the category below has finer brackets, to show variation on its own:
    # PublicInsuranceProportion_category_variation=case_when(
    #   PublicInsuranceProportion<0.1 ~"<0.1",
    #   PublicInsuranceProportion<0.2 ~"0.1-0.2",  
    #   PublicInsuranceProportion<0.3 ~"0.2-0.3",
    #   PublicInsuranceProportion<0.4 ~ "0.3-0.4",
    #   PublicInsuranceProportion<0.6 ~ "0.4-0.6",
    #   PublicInsuranceProportion<0.8 ~ "0.6-0.8",
    #   PublicInsuranceProportion<1.0 ~ "0.8-1.0",
    #   TRUE ~ 'No data'
    # )
  )
#map the data out, check area with higher public coverage:
# ggplot()+
#   geom_sf(data = year2019, aes( fill=PublicInsuranceProportion)) + 
#   scale_fill_viridis_c(direction=-1)+   # color shows from lightest to darkest, if 1(defualt), it's from darkest to lightest
#   xlab("Easting")+
#   ylab("Northing")+
#   ggtitle(paste0("Proportion of public insurance coverage by tract ","\n",
#                  "Reference: ACS 5 year 2015-2019"))

# #% public insurance coverage
# ggplot()+
#   geom_sf(data = year2019, aes( fill= PublicInsuranceProportion_category)) + 
#   # scale_fill_viridis_c(direction=-1)+   # color shows from lightest to darkest, if 1(defualt), it's from darkest to lightest
#   scale_fill_manual(values = c("<0.2" = "#7B886B", "0.2-0.4" = "#88BB92", "0.4-0.6" = "#94DDBC","0.6-0.8"="#A0ECD0","0.8-1.0"="white",
#                                "No data"="black"))+ # Assign specific colors to categories
#   
#   xlab("Easting")+
#   ylab("Northing")+
#   ggtitle(paste0("Proportion of public insurance coverage by tract ","\n",
#                  "Reference: ACS 5 year 2019"))
# #% DHIN representation
# ggplot()+
#   geom_sf(data = year2019, aes( fill= DHIN_representation_rate_category)) + 
#   # scale_fill_viridis_c(direction=-1)+   # color shows from lightest to darkest, if 1(defualt), it's from darkest to lightest
#   scale_fill_manual(values = c("<0.2" = "#7B886B", "0.2-0.4" = "#88BB92", "0.4-0.6" = "#94DDBC","0.6-0.8"="#A0ECD0","0.8-1.0"="white",
#                                "No data"="black"))+ # Assign specific colors to categories
#   
#   xlab("Easting")+
#   ylab("Northing")+
#   ggtitle(paste0("Proportion of DHIN representation by tract ","\n",
#                  "Reference: ACS 5 year 2019"))



final_data_before_simulation<-rbind(year2019,year2022)%>%mutate(
  
  shape1=FemaleOver40_wt_public_coverage+1   ,  #beta_i+1  , beta_i is the tract i's population who public insurance coverage (age>40 female), and is used for a beta distribution.
  shape2= CensusFemaleOver40-FemaleOver40_wt_public_coverage+1  #N_i - beta_i +1    , N_i is tract leve eligible women from census
)


#save datasets that are ready for adjustment
save.image(paste0(path,"Final_data_before_simulation_complete.RData"))
save(final_data_before_simulation,final_data_before_simulation_2019,final_data_before_simulation_2022,file=paste0(path,"Final_data_before_simulation.RData"))
