#Functions related to calculating RI from retention time indices

#1. Function to calculate RI from alkane RTs
alkane_RT_to_RI_calc<-function(	ms_data= raw_data,    #MSExperiment object from XCMS
                                alkanes= "/Users/diwalke/Dropbox/RESEARCH/Scripts/2024/3-Complete_GC-HRMS_processing_wf/10-Update_240608/3-Extraction_inputs/240806_Test_Batches_Alkanes.xlsx",
                                average_rt= TRUE,
                                filter_rt= TRUE,
                                mapfile= mapfile) {

  #Loop to adjust retention time with retention index
  #Results vector
  RI_rt_FINAL<- c() 	#Save adjusted RTs
  batch_n<- 0			#Batch column adjuster for alkanes_RT


  alkane_rt<-read_xlsx(alkanes)	#Read in alkane retention times

  #Remove any NA rows from alkane rt
  alkane_rt[alkane_rt == "NA"]<- NA
  alkane_rt<-alkane_rt[complete.cases(alkane_rt[, -c(1:5)]), ]

  alkane_rt<- alkane_rt %>%
    mutate_at(c(1, 6:ncol(alkane_rt)), ~ as.numeric(.))


  #If RTs average, first calculate average for each batch
  if(average_rt) {
    avg_rts<- c()
    for(kk in seq(6, ncol(alkane_rt), 2)) {
      avg_rts<- cbind(avg_rts, rowMeans(alkane_rt[, c(kk,kk+1)]))

    }
    alkane_rt<-cbind(alkane_rt[, 1:5], avg_rts)
  }

  #Determine min and max retention time
  min_time<- max(alkane_rt[1, -c(1:5)]) + 0.005
  max_time<- min(alkane_rt[nrow(alkane_rt), -c(1:5)]) - 0.005

  #Filter retention times to only include features with specified range
  if(filter_rt){
    ms_data<- filterRt(ms_data, rt = c(min_time * 60, max_time * 60))
  }

  #for loop to update retention times with retention index
  for(ii in unique(mapfile$Batch)) {

    #List of FileID numbers from sequence file
    batch_rownames<- as.numeric(rownames(mapfile[mapfile$Batch == ii, ]))

    #Select features detected in each file from corresponding batch
    raw_rt<- rtime(ms_data[batch_rownames])
    RI_rt_ii<- rep(0, length(raw_rt))

    #Loop through each pair of alkanes and calculate RI for each alkane pair
    for(jj in 2:nrow(alkane_rt)) {

      #Alkane retention times
      alkane_1<- as.numeric(alkane_rt[jj - 1, 6 + batch_n] * 60) #Tz
      alkane_2<- as.numeric(alkane_rt[jj, 6 + batch_n] * 60)    #Tz+1

      #Select all retention times between alkane_1 and alkane_2
      rt_jj<- which(raw_rt > alkane_1 & raw_rt <= alkane_2)

      #Calculate Retention time index
      RI_rt_ii[rt_jj]<- 100 * ((raw_rt[rt_jj] - alkane_1) / (alkane_2 - alkane_1) + as.numeric(alkane_rt[jj - 1, 1]))
    }
    batch_n<- batch_n + 1					#Advance to next batch RIs
    RI_rt_FINAL<-c(RI_rt_FINAL, RI_rt_ii)	#Final adjusted RTs
  }
  print("Retention time indices calculation complete")

  #Update raw retention time with alkane indices
  ms_data@spectra$rtime<- RI_rt_FINAL

  return(ms_data)

} #End of function 1



