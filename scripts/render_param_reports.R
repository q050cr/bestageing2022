

rm(list=ls())  # not in rmd doc, otherwise params deleted

## PARAMS -----------------------------------------------------------
# should modeling be performed again?
model.new <- FALSE

## REPORT and MAIL
library(dplyr)
gmailr::gm_auth_configure(path = "../library/credentials.json")
gmailr::gm_auth(email=TRUE, cache = ".secret")

render_report = function(disease, analysis) {
  ## disease: acs, cad, hfref, dcm | default: dcm
  ## analysis: selected, full | default: selected
  
  rmarkdown::render(
    "scripts/main.Rmd",
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
    output_dir = "reports",
    envir = globalenv()
  )
  # save filename for sending later
  filename.html <<-  paste0("reports/",
                            format(Sys.time(), "%Y%m%d"), '_ba_diagnostic_DISEASE_', stringr::str_to_upper(disease),
                            "_ANALYSIS_", stringr::str_to_upper(analysis), '.html'
  )
  title.mail <<- paste0("Diagnostic performance of miRNAs in ", stringr::str_to_upper(disease) )
}

# render scripts ----------------------------------------------------------
diseases <- c("dcm", "acs", "cad", "hfref")
analysis <- c("selected", "full")
all_combis <- tidyr::crossing(diseases, analysis) %>% 
  filter(analysis=="selected")# & diseases=="dcm")# %>% 

if (model.new == TRUE ) {
  source(file = "scripts/model.R")
}

time11 <- Sys.time()
# render for all diseases both selected and all miRNAs analysis
for (reports in 1:nrow(all_combis)) {
  render_report(disease = all_combis$diseases[reports], analysis = all_combis$analysis[reports])
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


