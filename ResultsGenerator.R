

filelist <- list.files("./data")

# loop through the file list to read in data and clean it up

for (results_file in filelist) {
  
  
  rmarkdown::render(input = "genericreport.rmd", 
                    output_format = 'pdf_document',
                    output_file = paste(results_file, ".pdf", sep=''),
                    output_dir = "reports",
                    params = list(
                      set_title = "WooTitle",
                      set_date = Sys.Date()
                    ))
  
}