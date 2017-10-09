library(readxl)
library(stringr)

directory <- "./data/"
filelist <- list.files(directory)
report_dir<- "reports"


element_list<-c("C","H","N","O","S","F","Cl","Br","I","Si")

element_count<-function(formula, element){
  if(is.na(formula)){
    return(0)
  }
  count<-str_match(formula, paste0(element,"([0-9]*)"))[2]
  if(is.na(count)){
    return(0)
  } else if (count==""){
    count<-"1"
  }
  return(as.numeric(count))
}

V_element_count<-Vectorize(element_count, vectorize.args = c("formula"), SIMPLIFY = TRUE)

summarize_table<-function(ecount, ct=compound_table){
  cut_table<-ct[ct$C == ecount[1] & ct$H == ecount[2] & ct$N==ecount[3] &
                  ct$O == ecount[4] & ct$S == ecount[5] & ct$F == ecount[6] &
                  ct$Cl == ecount[7] & ct$Br == ecount[8] & ct$I == ecount[9] &
                  ct$Si == ecount[10],]
  return(sum(as.numeric(cut_table$Area)))
}

# loop through the file list to read in data and clean it up

for (results_file in filelist) {
  
  results<-as.data.frame(read_excel(paste0(directory,results_file), col_names = FALSE))
  
  ctrow<-which(results[,1]=='Compound Table')
  
  header<-results[1:(ctrow-1),]
  compound_table<-results[(ctrow+1):nrow(results),]
  
  header<-Filter(function(x)!all(is.na(x)), header)
  compound_table <- Filter(function(x)!all(is.na(x)), compound_table)

  header_vals<-data.frame("key"=character(), "value"=character())
  
  for(i in seq(1, ncol(header), by = 2)){
    h<-unname(header[,i:(i+1)])
    h<-h[complete.cases(h),]
    colnames(h)<-c("key","value")
    header_vals<-rbind(header_vals, h)
  }
  
  rm(header, h, results)
  
  data_file_name<-header_vals[header_vals$key=="Data File", 2]
  sample_id<-header_vals[header_vals$key=="Sample Name", 2]
  sample_name<-header_vals[header_vals$key=="Comment", 2]
  analysis_date<-header_vals[header_vals$key=="Acquired Time", 2]
  acquisition_method<-header_vals[header_vals$key=="Acq Method", 2]
  data_analysis_method<-header_vals[header_vals$key=="DA Method", 2]
  operator_name<-header_vals[header_vals$key=="User Name", 2]
  
  colnames(compound_table)<-compound_table[1,]
  compound_table<-compound_table[2:(nrow(compound_table)-2),]
  
  for(i in 1:length(element_list)){
    compound_table[element_list[i]]<-V_element_count(formula = compound_table$`Molecular Formula`, element = element_list[i])
  }
  
  etable<-compound_table[,colnames(compound_table) %in% element_list]
  
  total_area<-sum(as.numeric(compound_table$Area))
  
  by_all <- unique(etable)
  rownames(by_all) <- NULL
  
  by_all$Area<-apply(by_all, 1, function(x) summarize_table(x, compound_table))
  by_all$`Area%` <- (by_all$Area/total_area)*100
  
  by_c<-data.frame("Carbon Count"=integer(), "Area"=integer(), "Area%"=numeric())
  by_ch<-data.frame("Carbon Count"=integer(), "Hydrogen Count"=integer(), "Area"=integer(), "Area%"=numeric())

  for(i in unique(compound_table$C)){
    ct<-compound_table[compound_table$C == i,]
    nr<-data.frame("Carbon Count"=i, "Area"=sum(as.numeric(ct$Area)), "Area%"=(sum(as.numeric(ct$Area))/total_area)*100)
    by_c<-rbind(by_c, nr)
    for(j in unique(ct$H)){
      ch<-ct[ct$H == j,]
      nh<-data.frame("Carbon Count"=i, "Hydrogen Count"=j, "Area"=sum(as.numeric(ch$Area)), "Area%"=(sum(as.numeric(ch$Area))/total_area)*100)
      by_ch<-rbind(by_ch, nh)
    }
  }
  
  by_c[order(by_c$Carbon.Count),]
  by_ch[order(by_ch$Carbon.Count, by_ch$Hydrogen.Count), ]
  by_all[order(by_all$C, by_all$H, by_all$N, by_all$O, by_all$S, by_all$F, by_all$Cl, by_all$Br, by_all$I, by_all$Si),]
  
  rmarkdown::render(input = "genericreport.rmd", 
                    output_format = 'pdf_document',
                    output_file = paste(sample_id, "_breakdown.pdf", sep=''),
                    output_dir = report_dir,
                    params = list(
                      set_title = "WooTitle",
                      set_date = Sys.Date()
                    ))
  
  write.csv(x = by_c, file = paste0(report_dir, "/", sample_id, "_by_c.csv"), row.names = FALSE)
  write.csv(x = by_ch, file = paste0(report_dir, "/", sample_id, "_by_ch.csv"), row.names = FALSE)
  write.csv(x = by_all, file = paste0(report_dir, "/", sample_id, "_by_all.csv"), row.names = FALSE)
  
}


