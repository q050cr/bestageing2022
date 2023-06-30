

# FNS for gpt api -----------------------------------------------------------

# CONSIDER
## Default rate limits  https://platform.openai.com/docs/guides/rate-limits
## https://github.com/openai/openai-cookbook/blob/main/examples/How_to_handle_rate_limits.ipynb
## 3,500 RPM (requests per minute)
## 90,000 TPM (tokens per minute)


library(openai)
library(chatgpt)
library(httr)
library(jsonlite)
library(reticulate)  # for 'tiktoken' library to count tokens
library(tidyverse)
library(glue)
use_python("/Library/Frameworks/Python.framework/Versions/3.11/bin/python3", required = T)
tiktoken <- import("tiktoken")
encoding <- tiktoken$encoding_for_model("gpt-3.5-turbo")  # model

#myapikey <- Sys.getenv("OPENAI_API_KEY")
#gptmodel <- openai::list_models()$data$root[39]   # gpt-3.5-turbo-0613


## calc api cost per request ----------------------------------
calculate_input_token_cost <- function(abstract) {
  # fn returns a list containing $input_tokens and $input_cost
  cost_per_1000_input_tokens <- 0.0015  # https://openai.com/pricing 2023-06-15 dollars
  input_tokens <- length(encoding$encode(abstract))
  
  # Calc
  input_cost <- (input_tokens / 1000) * cost_per_1000_input_tokens
  return(list(input_tokens=input_tokens, input_cost=input_cost))
}

## chatGPT function ;) ------------------------------

# Function to check if a segment is JSON formatted
is_json_formatted <- function(segment) {
  grepl("^\\{.*\\}$", segment)
}

# https://www.listendata.com/2023/05/chatgpt-in-r.html
chatGPT <- function(
    system_prompt, 
    user_prompt,
    #abstract_mirna_list,  # not needed any more in fn
    cost_per_1000_input_tokens = 0.0015,  # as of 2023-06-15
    modelName = "gpt-3.5-turbo",
    temperature = 0.7,  # default is 1
    max_tokens = 2048,
    top_p = 1,  # default
    apiKey = Sys.getenv("OPENAI_API_KEY")) {
  ## returns a list containing 
  ## 4 variables: `usage_total_tokens`, `query_cost`, `time`, and the plain `answer` of the request
  
  # Parameters for gpt
  params <- list(
    model = modelName,
    temperature = temperature,  # https://gptforwork.com/guides/openai-gpt3-temperature
    max_tokens = max_tokens,
    top_p = top_p
  )
  
  if(nchar(apiKey)<1) {
    apiKey <- readline("Paste your API key here: ")
    Sys.setenv(OPENAI_API_KEY = apiKey)
  }
  
  # messages
  messages = list(
    list(
      role = "system",
      content = system_prompt
    ),
    list(
      role = "user", 
      content = user_prompt
    )
  )
  # did not work in R, array error...
  functions = list(
    list(
      "name" = "get_current_weather", 
      "description" = "Get the current weather in a given location",
      "parameters" = list(
        "type" = "object",
        "properties" = list(
          "location" = list(
            "type" = "string",
            "description" = "The city and state, e.g. San Francisco, CA"
          ),
          "unit" = list(
            "type" = "string",
            "enum" = list(
              "celsius, fahrenheit"
            )
          )
        ),
        "required" = "location"
      )
    )
  )
  
  # https://medium.com/dev-bits/a-clear-guide-to-openai-function-calling-with-python-dcbc200c5d7
  
  response <- POST(
    url = "https://api.openai.com/v1/chat/completions", 
    add_headers(Authorization = paste("Bearer", apiKey)),
    content_type_json(),
    encode = "json",
    # body = toJSON(c(params, list(messages = messages)), auto_unbox = TRUE)
    body = c(params, 
             list(
               messages = messages
             )
    )
  )
  
  if (response$status_code > 200) {
    stop(content(response))
  }
  
  time <- response$times[[6]]  # total time in seconds
  response <- content(response)
  usage_total_tokens <- response$usage$total_tokens
  # pricing
  query_cost <- (usage_total_tokens / 1000) * cost_per_1000_input_tokens  # in dollars
  
  answer <- trimws(response$choices[[1]]$message$content)  # or cat(content(response)$choices[[1]]$message$content)
  # cat(answer)
  
  return(list(answer=answer, usage_total_tokens=usage_total_tokens, 
              query_cost=query_cost, time=time))
}


# RETRY IF ERROR -----------------------------------------------------------
# Function that retries the API call until it succeeds
# message = "That model is currently overloaded with other requests. You can retry your request, or contact us through our help center at help.openai.com if the error persists
retry_chatGPT <- function(system_prompt, user_prompt, retries=3) {
  tryCatch({
    resultsFromQuery <- chatGPT(system_prompt = system_prompt, user_prompt = user_prompt, cost_per_1000_input_tokens = 0.0015,
                                modelName = "gpt-3.5-turbo", temperature = 0.7)
    return(resultsFromQuery)
  }, error = function(e) {
    if(retries > 0) {
      message("Error: ", conditionMessage(e), ". Retrying in 20 seconds...")
      Sys.sleep(20) # Optional: wait before retrying
      resultsFromQuery <- retry_chatGPT(system_prompt, user_prompt, retries - 1)
      return(resultsFromQuery)
    } else {  # after trying X times
      message("Error: ", conditionMessage(e), ". No retries left.")
      return(NULL)
    }
  })
}


# CLEAN JSON -----------------------------------------------------------
clean_text2json <- function(text){
  
  # Initialize empty list to store data
  data_list <- list()
  
  #  split the text into segments based on the blank lines
  segments <- strsplit(text, "\n\n")[[1]]
  miR_segments <- subset(segments, grepl("\\b(miR|let)-[0-9]+[a-z]?(-[0-9]+p)?\\*?\\b", segments))
  miR_names <- unlist(regmatches(miR_segments, gregexpr("\\b(miR|let)-[0-9]+[a-z]?(-[0-9]+p)?\\*?\\b", miR_segments)))
  
  # Loop through the text elements
  for (i in seq_along(segments)) {
    # Extract miRNA name
    miRNA_name <- str_extract(segments[i], "\\b(miR|let)-[0-9]+[a-z]?(-[0-9]+p)?\\*?\\b")
    
    # Extract JSON strings
    json_strs <- str_extract_all(segments[i], "\\{[^\\}]*\\}")[[1]]
    
    # Parse JSON and add miRNA name
    json_data <- lapply(json_strs, function(x) c(miRNA = miRNA_name, fromJSON(x)))
    
    json_data[[1]]$miRNA
    json_data[[1]]$related_topic
    json_data[[1]]$direction_upreg_downreg
    
    
    # Append to list
    data_list <- c(data_list, json_data)
  }
  
  # Combine list elements into dataframe
  df <- bind_rows(data_list)
  
  # View the dataframe
  print(df)
  
}
  
#   tryCatch({
#     if (length(miR_names) == length(jsonsegments)) {
#       for (seg in jsonsegments) {
#         # extract the name (e.g., 'For miR-1')
#         #name <- gsub("\n.+", "", seg)
#         #name <- gsub("For ", "", name)
#         #name <- sub(":$", "", name)  # remove trailing colon
#         # extract the JSON part from the string
#         json_part <- gsub("^.+:\n", "", seg)
#         # convert JSON string to a list
#         json_list <- fromJSON(json_part, simplifyDataFrame = TRUE)
#         # convert the list to a data frame
#         json_df <- as.data.frame(t(json_list))
#         
#         # add the name column to the data frame
#         #json_df$Name <- name
#         
#         # bind the new data frame with the main data frame
#         df_temp <- rbind(df_temp, json_df)
#       }
#       # view the data frame
#       df_temp_abstract_x <- unnest(as_tibble(df_temp), 
#                                    cols = c(related_topic, direction_upreg_downreg, primary_literature, serum_plasma_tissue, mortality, measurement_type, sample_size))
#       df_temp_abstract_x$miR_names <- miR_names
#       
#       # return
#       return(list(df_temp_abstract_x=df_temp_abstract_x, usage_total_tokens=usage_total_tokens, 
#                   query_cost=query_cost, time=time, answer=answer))
#       
#     } else {
#       print(glue("Length of queried 'miR_names' = {length(miR_names)} is 
#                    not equal to the length of queried json 'jsonsegments' = {length(jsonsegments)}, 
#                    which indicates that GPT-3.5 output was not processed correctly.")
#       )
#       return(df_temp_abstract_x=tibble(), usage_total_tokens=usage_total_tokens, 
#              query_cost=query_cost, time=time, answer=answer)
#     }
#   }, error = function(e) {
#     # Error handling code
#     # This code will run if an error occurs within the try block
#     print(paste("An error occurred:", e$message))
#     return(df_temp_abstract_x=tibble(), usage_total_tokens=usage_total_tokens, 
#            query_cost=query_cost, time=time, answer=answer)
#   }, warning = function(w) {
#     # Warning handling code
#     # This code will run if a warning occurs within the try block
#     # Replace the comment with your desired warning handling code
#     print(paste("A warning occurred:", w$message))
#     return(df_temp_abstract_x=tibble(), usage_total_tokens=usage_total_tokens, 
#            query_cost=query_cost, time=time, answer=answer)
#   })
# }
#   
#   


# create new "clean-json-function" ------------------------------------------
#   
#   
#   ### clean output -------------------------------------------
#   # split the text into segments based on the blank lines
#   segments <- strsplit(answer, "\n\n")[[1]]
#   miR_segments <- subset(segments, grepl("\\b(miR|let)-[0-9]+[a-z]?(-[0-9]+p)?\\*?\\b", segments))
#   miR_names <- unlist(regmatches(miR_segments, gregexpr("\\b(miR|let)-[0-9]+[a-z]?(-[0-9]+p)?\\*?\\b", miR_segments)))
#   
#   # Initialize an empty list to store JSON data
#   data_list <- list()
#   # Loop through the text elements
#   for (i in seq_along(segments)) {
#     # Extract miRNA name
#     miRNA_name <- str_extract(segments[i], "\\b(miR|let)-[0-9]+[a-z]?(-[0-9]+p)?\\*?\\b")
#     
#     # Extract JSON strings !!
#     json_strs <- str_extract_all(segments[i], "\\{[^\\}]*\\}")[[1]]
#     
#     # Parse JSON and add miRNA name
#     json_data <- lapply(json_strs, function(x) c(miRNA = miRNA_name, fromJSON(x, simplifyDataFrame = TRUE)))
#     
#     # Append to list
#     data_list <- c(data_list, json_data)
#   }
#   
#   
#   # Extract JSON data using regular expression
#   json_data <- str_extract_all(segments, "\\{(?:[^{}]|(?R))*\\}")[[1]]
#   
#   # Remove segments that are not JSON formatted
#   jsonsegments <- segments[is_json_formatted(segments)]
#   df_temp <- data.frame()
#   
#   tryCatch({
#     if (length(miR_names) == length(jsonsegments)) {
#       for (seg in jsonsegments) {
#         # extract the name (e.g., 'For miR-1')
#         #name <- gsub("\n.+", "", seg)
#         #name <- gsub("For ", "", name)
#         #name <- sub(":$", "", name)  # remove trailing colon
#         # extract the JSON part from the string
#         json_part <- gsub("^.+:\n", "", seg)
#         # convert JSON string to a list
#         json_list <- fromJSON(json_part, simplifyDataFrame = TRUE)
#         # convert the list to a data frame
#         json_df <- as.data.frame(t(json_list))
#         
#         # add the name column to the data frame
#         #json_df$Name <- name
#         
#         # bind the new data frame with the main data frame
#         df_temp <- rbind(df_temp, json_df)
#       }
#       # view the data frame
#       df_temp_abstract_x <- unnest(as_tibble(df_temp), 
#                                    cols = c(related_topic, direction_upreg_downreg, primary_literature, serum_plasma_tissue, mortality, measurement_type, sample_size))
#       df_temp_abstract_x$miR_names <- miR_names
#       
#       # return
#       return(list(df_temp_abstract_x=df_temp_abstract_x, usage_total_tokens=usage_total_tokens, 
#                   query_cost=query_cost, time=time, answer=answer))
#       
#     } else {
#       print(glue("Length of queried 'miR_names' = {length(miR_names)} is 
#                    not equal to the length of queried json 'jsonsegments' = {length(jsonsegments)}, 
#                    which indicates that GPT-3.5 output was not processed correctly.")
#       )
#       return(df_temp_abstract_x=tibble(), usage_total_tokens=usage_total_tokens, 
#              query_cost=query_cost, time=time, answer=answer)
#     }
#   }, error = function(e) {
#     # Error handling code
#     # This code will run if an error occurs within the try block
#     print(paste("An error occurred:", e$message))
#     return(df_temp_abstract_x=tibble(), usage_total_tokens=usage_total_tokens, 
#            query_cost=query_cost, time=time, answer=answer)
#   }, warning = function(w) {
#     # Warning handling code
#     # This code will run if a warning occurs within the try block
#     # Replace the comment with your desired warning handling code
#     print(paste("A warning occurred:", w$message))
#     return(df_temp_abstract_x=tibble(), usage_total_tokens=usage_total_tokens, 
#            query_cost=query_cost, time=time, answer=answer)
#   })
# }
# 
