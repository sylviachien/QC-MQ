# app10 function.R

########## Function.1 ##########################################################################################
# x = pg, output = parsed_pg
process_pg <- function(x){
  pg = x
  bos_contam_list = c("CON__P00766", "CON__Q2KIG3", "CON__Q0VCM5", "CON__Q3SZ57", "CON__Q9N2I2", "CON__Q3SZH5", "CON__P28800", "CON__Q1A7A4", "CON__P41361", "CON__Q2YDI2", "CON__Q3Y5Z3", "CON__P81644", "CON__Q2KJ83", "CON__Q2KIT0", "CON__A2I7N3", "CON__Q3SZV7", "CON__Q2KJC7", "CON__Q3SZR3", "CON__Q28107", "CON__P02672", "CON__Q1RMN8", "CON__Q58D62", "CON__P06868", "CON__Q2KJF1", "CON__P02584", "CON__P02777", "CON__Q3SX14", "CON__P17697", "CON__Q6T181", "CON__P34955", "CON__P21752", "CON__Q32PJ2", "CON__Q28194", "CON__P00978", "CON__Q5XQN5", "CON__Q32PI4", "CON__Q9TTE1", "CON__Q2KIU3", "CON__P01044-1", "CON__P67983", "CON__Q28065", "CON__Q862S4", "CON__Q2KIF2", "CON__Q3SX28", "CON__Q0V8M9", "CON__Q148H6", "CON__Q29RQ1", "CON__Q95M17", "CON__P07224", "CON__Q2HJF0", "CON__Q2KIH2", "CON__A2I7N0", "CON__P12763", "CON__P17690", "CON__P02769", "CON__P02676", "CON__P50448", "CON__P01030", "CON__P01966", "CON__P00735", "CON__Q03247", "CON__Q3ZBS7", "CON__Q2UVX4", "CON__Q9TT36", "CON__Q28085", "CON__Q3SX09", "CON__P01045-1", "CON__Q3ZBD7", "CON__Q3MHN2", "CON__Q9TRI1", "CON__P15497", "CON__Q95121", "CON__Q05443", "CON__P02070", "CON__Q2KIS7", "CON__Q3MHH8", "CON__Q3T052", "CON__Q3KUS7", "CON__Q1RMK2", "CON__Q2TBQ1", "CON__Q05B55", "CON__A2I7N1", "CON__P04258", "CON__Q2KJ62", "CON__Q0IIK2", "CON__Q3MHN5", "CON__P02662", "CON__P02663", "CON__P02666", "CON__P02668", "CON__P31096", "CON__P02754", "CON__P00711", "CON__P62894", "CON__Q29443", "CON__ENSEMBL:ENSBTAP00000006074", "CON__ENSEMBL:ENSBTAP00000038329", "CON__REFSEQ:XP_001252647", "CON__ENSEMBL:ENSBTAP00000007350", "CON__ENSEMBL:ENSBTAP00000038253", "CON__ENSEMBL:ENSBTAP00000023402", "CON__ENSEMBL:ENSBTAP00000024466", "CON__ENSEMBL:ENSBTAP00000023055", "CON__ENSEMBL:ENSBTAP00000018229", "CON__ENSEMBL:ENSBTAP00000016046", "CON__ENSEMBL:ENSBTAP00000024462", "CON__ENSEMBL:ENSBTAP00000014147", "CON__ENSEMBL:ENSBTAP00000033053", "CON__ENSEMBL:ENSBTAP00000001528", "CON__ENSEMBL:ENSBTAP00000037665", "CON__ENSEMBL:ENSBTAP00000031900", "CON__ENSEMBL:ENSBTAP00000031360", "CON__ENSEMBL:ENSBTAP00000018574", "CON__ENSEMBL:ENSBTAP00000032840", "CON__ENSEMBL:ENSBTAP00000011227", "CON__ENSEMBL:ENSBTAP00000025008", "CON__ENSEMBL:ENSBTAP00000034412", "CON__ENSEMBL:ENSBTAP00000013050", "CON__ENSEMBL:ENSBTAP00000016285", "CON__ENSEMBL:ENSBTAP00000024146", "CON__REFSEQ:XP_585019")
  names = colnames(pg)
  # check the column of 'contaminant', whether its 'Contaminant' or 'Potential.contaminant'
  
  if(is.element('Potential.contaminant', names)){
    
    pg = pg %>% 
      select(Protein.IDs,
             Protein.names,
             Gene.names,
             Fasta.headers,
             Only.identified.by.site,
             Reverse,
             Potential.contaminant,
             starts_with("Intensity."),
             starts_with("LFQ."))
    parsed_pg = pg %>% 
      mutate(entry = strsplit(as.character(Protein.IDs), ";"))
    for(i in 1:nrow(parsed_pg)){
      parsed_pg$con_occurence[i] <- list(str_count(parsed_pg$entry[i],"CON_"))
    }
    parsed_pg <-  parsed_pg %>% 
      mutate(count_entry = str_count(Protein.IDs, ";")+1) %>% 
      mutate(uniq_contam = ifelse(count_entry == con_occurence, "+", ""))
    
    parsed_pg <- parsed_pg %>% 
      mutate(Bos.taurus = ifelse(is.element(entry, bos_contam_list), "+", ""))
    
    parsed_pg <- parsed_pg %>% 
      mutate(contam_kind = ifelse(Potential.contaminant!="+", "Non-Contaminant",
                                  ifelse(Potential.contaminant=="+" & Bos.taurus=="+" & uniq_contam=="+", "Contaminant.Bovine.Uniq",
                                         ifelse(Potential.contaminant=="+" & Bos.taurus=="+" & uniq_contam!="+", "Contaminant.Bovine.NonUniq",
                                                ifelse(Potential.contaminant=="+" & Bos.taurus!="+" & uniq_contam=="+", "Contaminant.NonBovine.Uniq",
                                                       ifelse(Potential.contaminant=="+" & Bos.taurus!="+" & uniq_contam!="+", "Contaminant.NonBovine.NonUniq", "NA"))))))
    parsed_pg <- parsed_pg %>% 
      mutate(kind_contam=ifelse(Protein.IDs=="CON__P00761","Trypsin",
                                ifelse(Protein.IDs=="CON__P02754", "LGB(Bovine)", contam_kind))) %>% 
      select(Protein.IDs, 
             Protein.names, 
             Gene.names, 
             Fasta.headers,
             Only.identified.by.site, 
             Reverse, 
             Potential.contaminant,
             Bos.taurus, 
             uniq_contam, 
             kind_contam, 
             count_entry,
             starts_with("Intensity."), 
             starts_with("LFQ."))
    
    return(parsed_pg)
    
  }
  if(is.element('Contaminant', names)){
    
    pg = pg %>% 
      select(Protein.IDs,
             Protein.names,
             Gene.names,
             Fasta.headers,
             Only.identified.by.site,
             Reverse,
             Contaminant,
             starts_with("Intensity."),
             starts_with("LFQ."))
    parsed_pg = pg %>% 
      mutate(entry = strsplit(as.character(Protein.IDs), ";"))
    for(i in 1:nrow(parsed_pg)){
      parsed_pg$con_occurence[i] <- list(str_count(parsed_pg$entry[i],"CON_"))
    }
    parsed_pg <-  parsed_pg %>% 
      mutate(count_entry = str_count(Protein.IDs, ";")+1) %>% 
      mutate(uniq_contam = ifelse(count_entry == con_occurence, "+", ""))
    
    parsed_pg <- parsed_pg %>% 
      mutate(Bos.taurus = ifelse(is.element(entry, bos_contam_list), "+", ""))
    
    parsed_pg <- parsed_pg %>% 
      mutate(contam_kind = ifelse(Contaminant!="+", "Non-Contaminant",
                                  ifelse(Contaminant=="+" & Bos.taurus=="+" & uniq_contam=="+", "Contaminant.Bovine.Unique",
                                         ifelse(Contaminant=="+" & Bos.taurus=="+" & uniq_contam!="+", "Contaminant.Bovine.NonUnique",
                                                ifelse(Contaminant=="+" & Bos.taurus!="+" & uniq_contam=="+", "Contaminant.NonBovine.Unique",
                                                       ifelse(Contaminant=="+" & Bos.taurus!="+" & uniq_contam!="+", "Contaminant.NonBovine.NonUnique", "NA"))))))
    parsed_pg <- parsed_pg %>% 
      mutate(kind_contam=ifelse(Protein.IDs=="CON__P00761","Trypsin",
                                ifelse(Protein.IDs=="CON__P02754", "LGB(Bovine)", contam_kind))) %>% 
      select(Protein.IDs, 
             Protein.names, 
             Gene.names, 
             Fasta.headers,
             Only.identified.by.site, 
             Reverse, 
             Contaminant,
             Bos.taurus, 
             uniq_contam, 
             kind_contam, 
             count_entry,
             starts_with("Intensity."), 
             starts_with("LFQ."))
    
    return(parsed_pg)
    
  }
}

#######################################################################################
#                                                                                     # 
#              Get Downloadable tables, input= x = parsed_pg, y = raw/lfq             #
#                                                                                     #  
#######################################################################################

########## Function.2 ##########################################################################################
getLog2Table <- function(x, y){
  names = colnames(x)
  if(y=="raw"){ extract = "Intensity."}
  if(y=="lfq"){ extract = "LFQ.intensity."}
  
  if(is.element('Potential.contaminant', names)){
    table = x %>% 
      filter(Reverse!="+") %>% 
      filter(Only.identified.by.site!="+") %>% 
      select(Protein.IDs, Protein.names, Gene.names, Potential.contaminant, kind_contam, starts_with(extract))
    n = ncol(table)
    table[,6:n] <- log(table[,6:n],2)
    table[table==-Inf] <- 0
    colnames(table) = gsub(extract,"",colnames(table))
    return(table)
  }
  if(is.element('Contaminant', names)){
    table = x %>% 
      filter(Reverse!="+") %>% 
      filter(Only.identified.by.site!="+") %>% 
      select(Protein.IDs, Protein.names, Gene.names, Contaminant, kind_contam, starts_with(extract))
    n = ncol(table)
    table[,6:n] <- log(table[,6:n],2)
    table[table==-Inf] <- 0
    colnames(table) = gsub(extract,"",colnames(table))
    return(table)
  }
}


########## Function.3 ##########################################################################################
getLog2CleanTable <- function(x, y){
  names = colnames(x)
  if(y=="raw"){ extract = "Intensity."}
  if(y=="lfq"){ extract = "LFQ.intensity."}
  
  if(is.element('Potential.contaminant', names)){
    table = x %>% 
      filter(Reverse!="+") %>% 
      filter(Only.identified.by.site!="+") %>% 
      filter(uniq_contam!="+") %>% 
      select(Protein.IDs, Protein.names, Gene.names, Potential.contaminant, kind_contam, starts_with(extract))
    n = ncol(table)
    table[,6:n] <- log(table[,6:n],2)
    table[table==-Inf] <- 0
    colnames(table) = gsub(extract,"",colnames(table))
    return(table)
  }
  if(is.element('Contaminant', names)){
    table = x %>% 
      filter(Reverse!="+") %>% 
      filter(Only.identified.by.site!="+") %>% 
      filter(uniq_contam!="+") %>% 
      select(Protein.IDs, Protein.names, Gene.names, Contaminant, kind_contam, starts_with(extract))
    n = ncol(table)
    table[,6:n] <- log(table[,6:n],2)
    table[table==-Inf] <- 0
    colnames(table) = gsub(extract,"",colnames(table))
    return(table)
  }
}


########## Function.4 ##########################################################################################
getLog2ContamTable <- function(x, y){
  names = colnames(x)
  if(y=="raw"){ extract = "Intensity."}
  if(y=="lfq"){ extract = "LFQ.intensity."}
  
  if(is.element('Potential.contaminant', names)){
    table = x %>% 
      filter(Potential.contaminant=="+") %>% 
      select(Protein.IDs, Protein.names, Gene.names, Potential.contaminant, kind_contam, starts_with(extract))
    n = ncol(table)
    table[,6:n] <- log(table[,6:n],2)
    table[table==-Inf] <- 0
    colnames(table) = gsub(extract,"",colnames(table))
    return(table)
  }
  if(is.element('Contaminant', names)){
    table = x %>% 
      filter(Contaminant=="+") %>% 
      select(Protein.IDs, Protein.names, Gene.names, Contaminant, kind_contam, starts_with(extract))
    n = ncol(table)
    table[,6:n] <- log(table[,6:n],2)
    table[table==-Inf] <- 0
    colnames(table) = gsub(extract,"",colnames(table))
    return(table)
  }
}


#######################################################################################
#                                                                                     # 
#                                Spike-in information                                 #
#                                                                                     #  
#######################################################################################

########## Function.5 ##########################################################################################
getSpikeInfo <- function(x, y){
  if(y=="raw"){ extract = "Intensity."}
  if(y=="lfq"){ extract = "LFQ.intensity."}
  #x[grepl("^P00000$", x$Protein.IDs),]
  s_data = x[grepl("^P00000$", x$Protein.IDs),] %>% 
  #s_data = x %>% filter(Protein.IDs=="P00000") %>%  
    select(starts_with(extract))
  data = t(s_data)
  data = data.frame(data, sample=rownames(data))
  colnames(data) <- c("intensity", "sample")
  data = data %>% mutate(log2.intensity = log(intensity,2))
  data$sample = gsub(extract,"", data$sample)
  data = data %>% select(sample, intensity, log2.intensity)
  return(data)
}

########## Function.6 ##########################################################################################
# x = parsed_pg, y = raw/lfq, z = protein.IDs
getCustomSpikeInfo <- function(x, y, z){
  if(y=="raw"){ extract = "Intensity."}
  if(y=="lfq"){ extract = "LFQ.intensity."}
  #x[grepl("^P00000$", x$Protein.IDs),]
  pname = paste0("^", z, "$")
  s_data = x[grepl(pname, x$Protein.IDs),] %>% 
    #s_data = x %>% filter(Protein.IDs=="P00000") %>%  
    select(starts_with(extract))
  data = t(s_data)
  data = data.frame(data, sample=rownames(data))
  colnames(data) <- c("intensity", "sample")
  data = data %>% mutate(log2.intensity = log(intensity,2))
  data$sample = gsub(extract,"", data$sample)
  data = data %>% select(sample, intensity, log2.intensity)
  return(data)
}

########## Function.7 ##########################################################################################
# get sample with spike below mean - 2sd, x=spike_dat, output:below
getSpikeBelow2sdSample <- function(x){
  m = mean(x$log2.intensity)
  s = sd(x$log2.intensity)
  below = x %>% filter(log2.intensity < (m-2*s))
  return(below)
}

########## Function.8 ##########################################################################################
# get sample with spike above mean + 2sd, x=spike_dat, output:above
getSpikeAbove2sdSample <- function(x){
  m = mean(x$log2.intensity)
  s = sd(x$log2.intensity)
  above = x %>% filter(log2.intensity > (m+2*s))
  return(above)
}

#######################################################################################
#                                                                                     # 
#                              contaminant information                                #
#                                                                                     #  
#######################################################################################

########## Function.9 ##########################################################################################
# get contaminant info, x= parsed_pg, output:raw/lfq_con_info
getContamInfo <- function(x, y){
  if(y=="raw"){ extract = "Intensity."}
  if(y=="lfq"){ extract = "LFQ.intensity."}
  
  info = x %>% 
    select(kind_contam, starts_with(extract)) %>% 
    gather(key="sample", value="intensity", -kind_contam)
  info$sample = gsub(extract, "", info$sample)
  return(info)
}

########## Function.10 ##########################################################################################
# get contaminant table, x= raw/lfq_con_info
getContamTable1 <- function(x){
  table = ddply(x, .(kind_contam, sample), summarise,
                sum_int=sum(intensity)) %>% 
    spread(kind_contam, sum_int)
  return(table)
}

########## Function.11 ##########################################################################################
# get contaminant table with sample info from doe, x= raw/lfq_con_info
getContamTable2 <- function(x, doe){
  info = x %>% 
    left_join(doe %>% select(sample.id, type, sample.type), by=c("sample"="sample.id"))
  table = ddply(info, .(kind_contam, sample.type), summarise,
                sum=sum(intensity),
                mean=mean(intensity),
                sd=sd(intensity),
                sem=sd(intensity)/sqrt(length(intensity))) %>%
    select(kind_contam,sample.type,mean) %>% 
    spread(kind_contam, mean)
  #return(table)
  return(info)
}


#######################################################################################
#                                                                                     # 
#                              protein count                                          #
#                                                                                     #  
#######################################################################################

########## Function.12 ##########################################################################################
# get protein count table for clean raw/lfq data, x=parsed_pg, output:prog_count
getProteinCount <- function(x, y){
  if(y=="raw"){ extract = "Intensity."}
  if(y=="lfq"){ extract = "LFQ.intensity."}
  
  table = x %>% 
    filter(Reverse!="+") %>% 
    filter(Only.identified.by.site!="+") %>% 
    filter(uniq_contam!="+") %>% 
    select(Protein.IDs, Protein.names, Gene.names, kind_contam, starts_with(extract))
  n = ncol(table)
  count = apply(table[,5:n], 2, function(x) sum(x>0))
  data = data.frame(count)
  data = data %>% 
    tibble::rownames_to_column("sample")
  data$sample = gsub(extract, "", data$sample)
  return(data)
}

########## Function.13 ##########################################################################################
# get protein count table from clean raw/lfq data with doe info, x=parsed_pg, y= raw/lfq
getProteinCountDoe <- function(x, y, doe){
  if(y=="raw"){ extract = "Intensity."}
  if(y=="lfq"){ extract = "LFQ.intensity."}
  
  table = x %>% 
    filter(Reverse!="+") %>% 
    filter(Only.identified.by.site!="+") %>% 
    filter(uniq_contam!="+") %>% 
    select(Protein.IDs, Protein.names, Gene.names, kind_contam, starts_with(extract))
  n = ncol(table)
  count = apply(table[,5:n], 2, function(x) sum(x>0))
  data = data.frame(count)
  data = data %>% 
    tibble::rownames_to_column("sample")
  data$sample = gsub(extract, "", data$sample)
  data = data %>% 
    left_join(doe, by=c("sample"="sample.id"))
  return(data)
}

########## Function.14 ##########################################################################################
# get expressed protein , x=parsed_pg, y=raw/lfq, output: per_data_melt_raw/lfq
getExpressedSamplePercentPerProteinEachType <- function(x, y, doe){
  if(y=="raw"){ extract = "Intensity."}
  if(y=="lfq"){ extract = "LFQ.intensity."}
  
  table = x %>% 
    filter(Reverse!="+") %>% 
    filter(Only.identified.by.site!="+") %>% 
    filter(uniq_contam!="+") %>% 
    select(Protein.IDs, Protein.names, Gene.names, kind_contam, starts_with(extract)) %>% 
    t()
  table <- data.frame(table[-c(1:4),])
  
  rownames(table) <- gsub(extract, "", rownames(table))
  
  table$sample.id <- rownames(table)
  table = table %>% 
    left_join(doe %>% select(sample.id, type))
  
  cnt <- table %>% filter(type=="control") %>% select(-type)
  cnt <- apply(cnt, 2, function(x) as.numeric(as.character(x)))
  per_cnt <- apply(t(cnt),1, function(x)(sum(x>0)*100/nrow(cnt)))
  
  smp <- table %>% filter(type=="sample") %>% select(-type)
  smp <- apply(smp, 2, function(x) as.numeric(as.character(x)))
  per_smp <- apply(t(smp),1, function(x)(sum(x>0)*100/nrow(smp)))
  
  per_data <- data.frame(Control=per_cnt, Sample=per_smp)
  per_data_melt <- melt(per_data)
  colnames(per_data_melt) <- c("Type", "Value")
  return(per_data_melt)
}

#######################################################################################
#                                                                                     # 
#                              peptide count                                          #
#                                                                                     #  
#######################################################################################

########## Function.15 ##########################################################################################
# get peptide count from pt x=pt, output:pep_count
getPepCount <- function(x){
  table = x %>% 
    select(starts_with("Intensity."))
  #  table <- apply(table, 2, function(x)as.numeric(x))
  count = apply(table, 2, function(x)sum(x>0))
  data = data.frame(count)
  data = data %>% 
    tibble::rownames_to_column("sample")
  data$sample = gsub("Intensity.", "", data$sample)
  return(data)
}

########## Function.16 ##########################################################################################
getPepCountDoe <- function(x, doe){
  table = x %>% 
    select(starts_with("Intensity."))
  #  table <- apply(table, 2, function(x)as.numeric(x))
  count = apply(table, 2, function(x)sum(x>0))
  data = data.frame(count)
  data = data %>% 
    tibble::rownames_to_column("sample")
  data$sample = gsub("Intensity.", "", data$sample)
  data = data %>% 
    left_join(doe, by=c("sample"="sample.id"))
  return(data)
}

#######################################################################################
#                                                                                     # 
#                              protein/peptide intensity                              #
#                                                                                     #  
#######################################################################################

########## Function.17 ##########################################################################################
# get total log2 protein intensity by individual sample x:table_log2_clean_raw/lfq
getTotalIntensity <- function(x){
  x <- x[,6:ncol(x)]
  Total.Log2.Intensity <- apply(x, 2, function(x)sum(x, na.rm = TRUE))
  data = data.frame(Total.Log2.Intensity)
  data = data %>% 
    tibble::rownames_to_column("sample")
  return(data)
}

########## Function.18 ##########################################################################################
# get log2 protein intensities by sample,  x:table_log2_clean_raw/lfq
getIntensity <- function(x){
  x <- x[,6:ncol(x)]
  data = data.frame(t(x))
  data = data %>% 
    tibble::rownames_to_column("sample")
  return(data)
}

########## Function.19 ##########################################################################################
# get total log2 peptide intensity by individual sample, x:pt output: pep_int_sum
getTotalIntensityPep <- function(x){
  x = x %>% 
    select(starts_with("Intensity."))
  n = ncol(x)
  x[,1:n] <- log(x[,1:n],2)
  x[x==-Inf] <- 0
  colnames(x) <- gsub("Intensity.", "", colnames(x))
  Total.Log2.Intensity <- apply(x, 2, function(x)sum(x, na.rm = TRUE))
  data = data.frame(Total.Log2.Intensity)
  data = data %>% 
    tibble::rownames_to_column("sample")
  return(data)
}

########## Function.20 ##########################################################################################
# median intensity plot:use pg_ints_raw/lfq_melt=x (table_log2_clean_r/l -> getIntensity -> melt), output:mp_raw/lfq
getMedianIntPlotData <- function(x, doe){
  doe <- doe %>% select(sample.id, type)
  x = x %>% 
    left_join(doe, by=c("sample"="sample.id"))
  dat_cal = x %>% 
    ddply(.(variable), summarise, sum=sum(value), mean=mean(value), sd=sd(value), sem=sd(value)/sqrt(length(value)))
  dat_mean_sum = sum(dat_cal$mean)
  dat_cal = dat_cal %>% 
    mutate(Proportion.Mean = mean/dat_mean_sum+0.000001)
  cnt = x %>% 
    filter(type=="control")
  cnt_cal = cnt %>% 
    ddply(.(variable), summarise, sum=sum(value), mean=mean(value), sd=sd(value), sem=sd(value)/sqrt(length(value)))
  cnt_mean_sum = sum(cnt_cal$mean)
  cnt_cal$mean = as.numeric(cnt_cal$mean)
  cnt_cal = cnt_cal %>% 
    mutate(Proportion.Mean = mean/cnt_mean_sum+0.000001)
  smp = x %>% 
    filter(type=="sample")
  smp_cal = smp %>% 
    ddply(.(variable), summarise, sum=sum(value), mean=mean(value), sd=sd(value), sem=sd(value)/sqrt(length(value)))
  smp_mean_sum = sum(smp_cal$mean)
  smp_cal$mean = as.numeric(smp_cal$mean)
  smp_cal = smp_cal %>% 
    mutate(Proportion.Mean = mean/smp_mean_sum+0.000001)
  merge = data.frame(#protein=dat_cal$variable,
    All=sort(as.numeric(dat_cal$Proportion.Mean)),
    Control=sort(as.numeric(cnt_cal$Proportion.Mean)),
    Sample=sort(as.numeric(smp_cal$Proportion.Mean)))
  return(merge)
}

#######################################################################################
#                                                                                     # 
#                                       top 20                                        #
#                                                                                     #  
#######################################################################################

########## Function.21 ##########################################################################################
# x=table_log2_clean_raw, output=melt_table_name_raw/lfq
getMeltwithNames <- function(x){
  names = colnames(x)
  if(is.element('Potential.contaminant', names)){
    data = x %>% 
      select(-Potential.contaminant, -kind_contam) %>% 
      gather(key="sample", value="value", -Protein.IDs, -Protein.names, -Gene.names)
    return(data)
  }
  if(is.element('Contaminant', names)){
    data = x %>% 
      select(-Contaminant, -kind_contam) %>% 
      gather(key="sample", value="value", -Protein.IDs, -Protein.names, -Gene.names)
    return(data)
  }
}

########## Function.22 ##########################################################################################
# x=melt_table_name_raw/lfq, output=annotated_mean_table_raw/lfq
countMeanbyProtein <- function(x){
  data = x %>% 
    ddply(.(Protein.IDs), summarise, Mean.Log2.Intensity=mean(value)) %>% 
    left_join(x %>% select(Protein.IDs, Protein.names, Gene.names), by="Protein.IDs") %>% 
    distinct() %>% 
    select(Protein.IDs, Protein.names, Gene.names, Mean.Log2.Intensity) %>% 
    arrange(desc(Mean.Log2.Intensity))
  data$Mean.Log2.Intensity = format(round(data$Mean.Log2.Intensity,2))
  data$Protein.IDs = gsub(";","; ",data$Protein.IDs)
  data$Protein.names = gsub(";","; ",data$Protein.names)
  data$Gene.names = gsub(";","; ",data$Gene.names)
  return(data)
}

########## Function.23 ##########################################################################################
## top20 sample.type
# x=melt_table_name_raw/lfq_doe
countMeanbyProteinSampletypeDoe <- function(x){
  info = x %>% select(Protein.IDs, Protein.names, Gene.names)
  data = x %>% 
    ddply(.(Protein.IDs, sample.type),  summarise, Mean.intensity=round(mean(value),2)) %>% 
    spread(key="sample.type", value="Mean.intensity")
  
  table = info %>% right_join(data) %>% 
    distinct()

  table$Protein.IDs = gsub(";","; ", table$Protein.IDs)
  table$Protein.names = gsub(";","; ", table$Protein.names)
  table$Gene.names = gsub(";", "; ", table$Gene.names)
  
  return(table)
}

########## Function.25 ##########################################################################################
## top20 type
# x=melt_table_name_raw/lfq_doe
countMeanbyProteinTypeDoe <- function(x){
  info = x %>% select(Protein.IDs, Protein.names, Gene.names)
  data = x %>% 
    ddply(.(Protein.IDs, type),  summarise, Mean.intensity=round(mean(value),2)) %>% 
    spread(key="type", value="Mean.intensity")
  
  table = info %>% right_join(data) %>% 
    distinct()
  
  table$Protein.IDs = gsub(";","; ", table$Protein.IDs)
  table$Protein.names = gsub(";","; ", table$Protein.names)
  table$Gene.names = gsub(";", "; ", table$Gene.names)
  
  return(table)
}




########## Function.23 ##########################################################################################
## top20 sample.type
# x=melt_table_name_raw/lfq_doe output:wide_sampletype_sort_table_raw/lfq
#countMeanbyProteinSampletypeDoe <- function(x){
#  data = x %>% 
#    ddply(.(Protein.IDs, sample.type),  summarise, Mean.intensity=mean(value)) %>% 
#    spread(key="sample.type", value="Mean.intensity")
#  # the columns are `Protein.IDs` and list of sample.type
#  # create separate tables for different sample.type
#  sampletype_name = colnames(data)[-1]
#  number_sampletype = ncol(data)-1
#  name_col = data %>% 
#    select(Protein.IDs) %>% 
#    left_join(x %>% select(Protein.IDs, Protein.names, Gene.names)) %>% 
#    distinct()
#  name_col$Protein.IDs = gsub(";","; ", name_col$Protein.IDs)
#  name_col$Protein.names = gsub(";","; ", name_col$Protein.names)
#  name_col$Gene.names = gsub(";", "; ", name_col$Gene.names)
#  big_table = ""
#  for(i in 1:number_sampletype){
#    j=i+1
#    table = data.frame(name_col, data[,j])
#    sort_table = table[order(-table[,4]),]
#    big_table = data.frame(big_table, sort_table)
#  }
#  big_table <- big_table[,-1]
#  for(i in 1:number_sampletype){
#    n = i*4
#    colnames(big_table)[n] <- sampletype_name[i]
#  } 
#  return(big_table)
#}

########## Function.24 ##########################################################################################
# x= wide_sampletype_sort_table_raw/lfq
getTop20sampletypeTable <- function(x){
  number_sampletype = ncol(x)/4
  sampletype_list = ""
  for(i in 1:number_sampletype){
    n=i*4
    sampletype_list = c(sampletype_list, colnames(x)[n])
  }
  sampletype_list = sampletype_list[-1]
  long_table = ""
  empty_col = ""
  m_x <- cbind(empty_col, x)
  for(i in 1:number_sampletype){
    n=i*4+1
    st_n = n-4
    small_table = m_x[1:20, st_n:n]
    small_table[,5] <- format(round(as.numeric(small_table[,5]),2))
    colnames(small_table) <- c("extra", "Protein.IDs", "Protein.names", "Gene.names", "Log2.Intensity")
    title_row = c("", paste0("Top 20 in sample.type: ",sampletype_list[i]), "", "", "")
    
    long_table = rbind(long_table, title_row, small_table)
  }
  long_table <- long_table[,-1]
  
  return(long_table)
}

########## Function.25 ##########################################################################################
## top20 type
# x=melt_table_name_raw/lfq_doe output:wide_sampletype_sort_table_raw/lfq
#countMeanbyProteinTypeDoe <- function(x){
#  data = x %>% 
#    ddply(.(Protein.IDs, type),  summarise, Mean.intensity=mean(value)) %>% 
#    spread(key="type", value="Mean.intensity")
#  # the columns are `Protein.IDs` and list of sample.type
#  # create separate tables for different sample.type
#  type_name = colnames(data)[-1]
#  number_type = ncol(data)-1
#  name_col = data %>% 
#    select(Protein.IDs) %>% 
#    left_join(x %>% select(Protein.IDs, Protein.names, Gene.names)) %>% 
#    distinct()
#  name_col$Protein.IDs = gsub(";","; ", name_col$Protein.IDs)
#  name_col$Protein.names = gsub(";","; ", name_col$Protein.names)
#  name_col$Gene.names = gsub(";", "; ", name_col$Gene.names)
#  big_table = ""
#  for(i in 1:number_type){
#    j=i+1
#    table = data.frame(name_col, data[,j])
#    sort_table = table[order(-table[,4]),]
#    big_table = data.frame(big_table, sort_table)
#  }
#  big_table <- big_table[,-1]
#  for(i in 1:number_type){
#    n = i*4
#    colnames(big_table)[n] <- type_name[i]
#  } 
#  return(big_table)
#}

########## Function.26 ##########################################################################################
# x= wide_sampletype_sort_table_raw/lfq
getTop20TypeTable <- function(x){
  number_type = ncol(x)/4
  type_list = ""
  for(i in 1:number_type){
    n=i*4
    type_list = c(type_list, colnames(x)[n])
  }
  type_list = type_list[-1]
  long_table = ""
  empty_col = ""
  m_x <- cbind(empty_col, x)
  for(i in 1:number_type){
    n=i*4+1
    st_n = n-4
    small_table = m_x[1:20, st_n:n]
    small_table[,5] <- format(round(as.numeric(small_table[,5]),2))
    colnames(small_table) <- c("extra", "Protein.IDs", "Protein.names", "Gene.names", "Log2.Intensity")
    title_row = c("", paste0("Top 20 in type: ",type_list[i]), "", "", "")
    
    long_table = rbind(long_table, title_row, small_table)
  }
  long_table <- long_table[,-1]
  
  return(long_table)
}

#######################################################################################
#                                                                                     # 
#                                       pca                                           #
#                                                                                     #  
#######################################################################################

########## Function.27 ##########################################################################################
# x = table_log2_clean_raw/lfq, output: pca_table_raw
getPCAtable <- function(x){
  x = data.frame(t(x))
  x = x[-c(1:5),]
  n = ncol(x)
  x[,1:n] <- lapply(x[,1:n],as.numeric)
  pca = prcomp(x)
  y = data.frame(pca$x[,1:2])
  y$sample = rownames(y)
  return(y)
}

########## Function.28 ##########################################################################################
# get contribution of PC1 and PC2
# x = table_log2_clean_raw/lfq, output: pc12_table_raw, a vector with 2 elements
getPC12Percentage <- function(x){
  x = data.frame(t(x))
  x = x[-c(1:5),]
  n = ncol(x)
  x[,1:n] <- lapply(x[,1:n],as.numeric)
  pca = prcomp(x)
  pc1_pv = round(unlist(pca$sdev^2/sum(pca$sdev^2))[1]*100, digits = 2)
  pc2_pv <- round(unlist(pca$sdev^2/sum(pca$sdev^2))[2]*100, digits = 2)
  y = c(pc1_pv, pc2_pv)
  return(y)
}

########## Function.29 ##########################################################################################
# x = table_log2_clean_raw/lfq, doe
getPCAtableDOE <- function(x, doe){
  x = data.frame(t(x))
  x = x[-c(1:5),]
  n = ncol(x)
  x[,1:n] <- lapply(x[,1:n],as.numeric)
  pca = prcomp(x)
  y = data.frame(pca$x[,1:2])
  y$sample = rownames(y)
  y = y %>% 
    left_join(doe %>% select(sample.id, type, sample.type, run.order), by=c("sample"="sample.id"))
  return(y)
}

#######################################################################################
#                                                                                     # 
#                            contribution of the variance                             #
#                                                                                     #  
#######################################################################################

########## Function.30 ##########################################################################################
# input x= table_log2_clean_raw/lfq, output: pvca_raw/lfq
getPVCAresult <- function(x, doe){
  
  d = doe
  rownames(d) = d$sample.id
  pdata = new("AnnotatedDataFrame", data=d)
  
  x = x[,-c(1:5)]
  x[x==-Inf] <- 0
  n = ncol(x)
  x[,1:n] <- lapply(x[,1:n],as.numeric)
  exprs = as.matrix(x)
  
  colnames(exprs) <- rownames(pdata)
  set = ExpressionSet(assayData = exprs,
                      phenoData = pdata)
  batch.factors = c("sample.type", "type")
  pvca.obj = pvcaBatchAssess(set, batch.factors, 0.6)
  
  result = data.frame(Variable = pvca.obj$label,
                      Proportion.Variance = pvca.obj$dat[1:4])
  return(result)
}

#######################################################################################
#                                                                                     # 
#                                   quantro test                                     #
#                                                                                     #  
#######################################################################################

########## Function.31 ##########################################################################################
# input x = table_log2_clean_raw/lfq, doe, output= qt_raw/lfq
doQuantro <- function(x, doe){
  x = x[,-c(1:5)]
  n = ncol(x)
  x[,1:n] <- lapply(x[,1:n],as.numeric)
  colnames(doe) <- tolower(colnames(doe))
  pdoe = doe %>% 
    arrange(sample.id)
  
  if(length(unique(pdoe$sample.type))<2){
    no_quantro <- "The number of 'sample.type' is less than 2, Quantro can not be applied!"
    return(no_quantro)
  }else{
    d = quantro(x, groupFactor=pdoe$sample.type, useMedianNormalized=FALSE, B=1000) 
    
    t = data.frame('Quantro_result' = c("nGroups", "nTotalSample", "quantroStat", "quantroPvaluePermutation"), 
                   'value'= c(unname(unlist(summary(d)[1])),
                              unname(unlist(summary(d)[2])),
                              round(quantroStat(d),digits=4),
                              ifelse(quantroPvalPerm(d)==0,"<0.001",round(quantroPvalPerm(d),digits=4))))
    return(t)
  }
  
}


########## Function.32 ##########################################################################################
# get anovar table from quantro result, x= table_log2_raw/lfq, output= qt_raw/lfq_a
doQuantroAnova <- function(x, doe){
  x = x[,-c(1:5)]
  x[x==-Inf] <- 0
  n = ncol(x)
  x[,1:n] <- lapply(x[,1:n],as.numeric)
  colnames(doe) <- tolower(colnames(doe))
  pdoe = doe %>% 
    arrange(sample.id)
  
  if(length(unique(pdoe$sample.type))<2){
    no_quantro <- "The number of 'sample.type' is less than 2, Quantro can not be applied!"
    return(no_quantro)
  }else{
    d = quantro(x, groupFactor=pdoe$sample.type, useMedianNormalized=FALSE, B=1000) 
    a = anova(d)
    return(a)
  }
  
}


########## Function.33 ##########################################################################################
# input x = table_log2_clean_raw/lfq, doe, output= qt_raw/lfq
doQuantroNoControl <- function(x, doe){
  x = x[,-c(1:5)]
  n = ncol(x)
  x[,1:n] <- lapply(x[,1:n],as.numeric)
  colnames(doe) <- tolower(colnames(doe))
  doe$type <- tolower(doe$type)
  
  pdoe = doe %>% 
    filter(type=="sample") %>% 
    arrange(sample.id)
  pdoe$sample.type <- factor(pdoe$sample.type)
  
  if(length(unique(pdoe$sample.type))<2){
    no_quantro <- "The number of 'sample.type' is less than 2, Quantro can not be applied!"
    return(no_quantro)
  }else{
    tt = as.data.frame(t(x))
    tt$sample.id = rownames(tt)
    m = tt %>% left_join(doe %>% select(sample.id, type)) %>% 
      filter(type=="sample") %>% 
      select(-type)
    xx = as.data.frame(t(m))
    colnames(xx) <- unlist(xx["sample.id",])
    xxx <- as.data.frame(xx[-nrow(xx),])
    j = ncol(xxx)
    xxx[,1:j] <- lapply(xxx[,1:j],as.character)
    xxx[,1:j] <- lapply(xxx[,1:j],as.double)
    
    d = quantro(xxx, groupFactor=pdoe$sample.type, useMedianNormalized=FALSE, B=1000) 
    
    t = data.frame('Quantro_result' = c("nGroups", "nTotalSample", "quantroStat", "quantroPvaluePermutation"), 
                   'value'= c(unname(unlist(summary(d)[1])),
                              unname(unlist(summary(d)[2])),
                              round(quantroStat(d),digits=4),
                              ifelse(quantroPvalPerm(d)==0,"<0.001",round(quantroPvalPerm(d),digits=4))))
    return(t)
  }
  
}


########## Function.34 ##########################################################################################
# get anovar table from quantro result, x= table_log2_raw/lfq, output= qt_raw/lfq_a
doQuantroNoControlAnova <- function(x, doe){
  x = x[,-c(1:5)]
  n = ncol(x)
  x[,1:n] <- lapply(x[,1:n],as.numeric)
  colnames(doe) <- tolower(colnames(doe))
  doe$type <- tolower(doe$type)
  
  pdoe = doe %>% 
    filter(type=="sample") %>% 
    arrange(sample.id)
  pdoe$sample.type <- factor(pdoe$sample.type)
  
  if(length(unique(pdoe$sample.type))<2){
    no_quantro <- "The number of 'sample.type' is less than 2, Quantro can not be applied!"
    return(no_quantro)
  }else{
    tt = as.data.frame(t(x))
    tt$sample.id = rownames(tt)
    m = tt %>% left_join(doe %>% select(sample.id, type)) %>% 
      filter(type=="sample") %>% 
      select(-type)
    xx = as.data.frame(t(m))
    colnames(xx) <- unlist(xx["sample.id",])
    xxx <- as.data.frame(xx[-nrow(xx),])
    j = ncol(xxx)
    xxx[,1:j] <- lapply(xxx[,1:j],as.character)
    xxx[,1:j] <- lapply(xxx[,1:j],as.double)
    
    d = quantro(xxx, groupFactor=pdoe$sample.type, useMedianNormalized=FALSE, B=1000) 
    
    t = data.frame('Quantro_result' = c("nGroups", "nTotalSample", "quantroStat", "quantroPvaluePermutation"), 
                   'value'= c(unname(unlist(summary(d)[1])),
                              unname(unlist(summary(d)[2])),
                              round(quantroStat(d),digits=4),
                              ifelse(quantroPvalPerm(d)==0,"<0.001",round(quantroPvalPerm(d),digits=4))))
    a = anova(d)
    return(a)
  }
  
  

}



# end

