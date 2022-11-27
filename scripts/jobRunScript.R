
# jobRunScript() ----------------------------------------------------------

rstudioapi::jobRunScript(path = "~/bestageing2022/scripts/render_param_reports.R", 
                         name = "jobscript_mirna_diagnostic",
                         workingDir="~/bestageing2022",
                         importEnv=FALSE,
                         exportEnv=""  # default (skip export)
                         )


# callr -------------------------------------------------------------------

callr::rscript( # https://callr.r-lib.org/reference/rscript.html
  script = "scripts/render_param_reports.R",
  cmdargs = character(),
  libpath = .libPaths()[1],
  #repos = default_repos(),
  #stdout = NULL,  # Optionally a file name to send the standard output to
  #stderr = NULL,
  #poll_connection = TRUE,
  echo = TRUE,  # default=FALSE
  show = TRUE,
  #callback = NULL,
  #block_callback = NULL,
  #spinner = FALSE,
  #system_profile = FALSE,
  #user_profile = "project",
  #env = rcmd_safe_env(),
  timeout = as.difftime(12, units = "hours"),  # default=Inf
  wd = "~/bestageing2022",  # defaults to the current working directory.
  fail_on_status = TRUE,  # default=FALSE
  color = TRUE,
  #...
)
