

## render report and mail results
library(gmailr)
library(dplyr)

rm(list=ls())  # not in rmd doc, otherwise params deleted

# define render function to stratify groups
render_report = function(disease, analysis) {
  ## disease: acs, cad, hfref, dcm | default: dcm
  ## analysis: selected, full | default: selected
  
  rmarkdown::render(
    "bestageing2022/scripts/20221110main.Rmd", 
    params = list(
      disease = disease,
      analysis = analysis,
      doc_title = paste0("Diagnostic performance of miRNAs in ", stringr::str_to_upper(disease) )
    ),
    output_file <<- paste0(
      format(Sys.time(), "%Y%m%d"), '_ba_diagnostic_DISEASE_', stringr::str_to_upper(disease),
      "_ANALYSIS_", stringr::str_to_upper(analysis), '.html'
      ),
    output_dir = "bestageing2022/reports",
    envir = globalenv()
  )
}

render_report(disease = "dcm", analysis = "selected")



diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>% 
  filter(analysis=="selected") %>% 
  slice(1:2)

## render first report
#time1 <- Sys.time()
#render_report(disease = "dcm", analysis = "selected")
#time2 <- Sys.time()
#time.diff.grid <- time2-time1

time11 <- Sys.time()
# render for all diseases both selected and all miRNAs analysis
for (i in 1:nrow(all_combis)) {
  render_report(disease = all_combis$diseases[i], analysis = all_combis$analysis[i])
  # send report
  email <-
   gm_mime() %>%
   gm_to("reich.c@hotmail.com") %>%
   gm_from("q050cr@gmail.com") %>%
   gm_subject(paste0("Job finished: ", output_file)) %>%
   gm_text_body("Hi Christoph,\nHere is a new report for you.\nCheers")
  email <- gm_attach_file(test_email, paste0("bestageing2022/reports/", output_file))
  gm_send_message(email)
}
time22 <- Sys.time()
time.diff.report <- time22-time11

time.diff.report

