library(readxl)
library(stringr)

filelist <- list.files("./data")


element_list<-c("C","H","N","O","S","F","Cl","Br","I","Si")

# loop through the file list to read in data and clean it up

for (results_file in filelist) {
  
  results<-as.data.frame(read_excel(results_file, col_names = FALSE))
  
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
  etable<-compound_table[,colnames(compound_table) %in% element_list]
  
  total_area<-sum(compound_table$area)
  
  for(i in 1:length(element_list)){
    compound_table[element_list[i]]<-V_element_count(formula = compound_table$`Molecular Formula`, element = element_list[i])
  }
  
  by_c<-data.frame("Carbon Count"=integer(), "Area"=integer(), "Area %"=numeric())
  by_ch<-data.frame("Carbon Count"=integer(), "Hydrogen Count"=integer(), "Area"=integer(), "Area %"=numeric())

  for(i in unique(compound_table$C)){
    ct<-compound_table[compound_table$C == i,]
    nr<-data.frame("Carbon Count"=i, "Area"=sum(ct$Area), "Area %"=(sum(ct$Area)/total_area)*100)
    by_c<-rbind(by_c, cf)
    for(j in unique(ct$H)){
      ch<-ct[ct$H == j,]
      nr<-data.frame("Carbon Count"=i, "Hydrogen Count"=j, "Area"=sum(ch$Area), "Area %"=(sum(ch$Area)/total_area)*100)
    }
  }
  
  by_all<- etable
  by_all$Area<-sapply(etable, function(x) element_count(x, compound_table))
  by_all$`Area %` <- by_all$Area/total_sum
  
  
  rmarkdown::render(input = "genericreport.rmd", 
                    output_format = 'pdf_document',
                    output_file = paste(sample_id, "breakdown.pdf", sep=''),
                    output_dir = "reports",
                    params = list(
                      set_title = "WooTitle",
                      set_date = Sys.Date()
                    ))
  
}

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
  return(as.integer(count))
}

V_element_count<-Vectorize(element_count, vectorize.args = c("formula"), SIMPLIFY = TRUE)

summarize_table<-function(ecount, ct=compound_table){
  cut_table<-ct[ct$C == ecount$C & ct$H == ecount$H & ct$N==ecount$N &
                ct$O == ecount$O & ct$S == ecount$S & ct$F == ecount$F &
                ct$Cl == ecount$Cl & ct$Br == ecount$Br & ct$I == ecount$I &
                ct$Si == ecount$Si,]
  return(sum(cut_table$Area))
}
