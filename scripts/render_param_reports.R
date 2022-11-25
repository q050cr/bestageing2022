



rm(list=ls())  # not in rmd doc, otherwise params deleted

# change libpath
myPaths <- .libPaths()
myPaths <- c("library", myPaths)[1:2]  # only use ow
.libPaths(myPaths)

## render report and mail results
library(dplyr)
gmailr::gm_auth_configure(path = "library/credentials.json")
gmailr::gm_auth(email=TRUE, cache = ".secret")

render_report = function(disease, analysis) {
  ## disease: acs, cad, hfref, dcm | default: dcm
  ## analysis: selected, full | default: selected
  
  rmarkdown::render(
    "bestageing2022/scripts/20221110main.Rmd",
    params = list(
      disease = disease,
      analysis = analysis,
      doc_title =  paste0("Diagnostic performance of miRNAs in ", stringr::str_to_upper(disease) )
    ),
    output_file = paste0(
      format(Sys.time(), "%Y%m%d"), '_ba_diagnostic_DISEASE_', stringr::str_to_upper(disease), 
      "_ANALYSIS_", stringr::str_to_upper(analysis),
      '.html'
    ),
    output_dir = "bestageing2022/reports",
    envir = globalenv()
  )
  # save filename for sending later
  filename.html <<-  paste0("bestageing2022/reports/",
                            format(Sys.time(), "%Y%m%d"), '_ba_diagnostic_DISEASE_', stringr::str_to_upper(disease),
                            "_ANALYSIS_", stringr::str_to_upper(analysis), '.html'
  )
  title.mail <<- paste0("Diagnostic performance of miRNAs in ", stringr::str_to_upper(disease) )
}

# render_report(disease = "dcm", analysis = "selected")
# render scripts ----------------------------------------------------------

diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>% 
  filter(analysis=="full" & diseases=="dcm")# %>% 
  #slice(1:2)

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
   gmailr::gm_mime() %>%
   gmailr::gm_to("reich.c@hotmail.com") %>%
   gmailr::gm_from("q050cr@gmail.com") %>%
   gmailr::gm_subject(paste0("Job finished: ", title.mail)) %>%
   gmailr::gm_text_body("Hi Christoph,\nHere is a new report for you.\nCheers")
  # attachment
  email <- gmailr::gm_attach_file(email, filename.html)
  gmailr::gm_send_message(email, user_id="me")
}
time22 <- Sys.time()
time.diff.report <- time22-time11

time.diff.report

