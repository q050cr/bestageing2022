
import os
import json
import pickle
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)  # Display any number of columns
pd.set_option('display.expand_frame_repr', False)  # Do not wrap to multiple pages

import openai
import time
from datetime import datetime
# import requests
# from tenacity import retry, wait_random_exponential, stop_after_attempt

# set variables
GPT_MODEL = "gpt-3.5-turbo-0613"
openai.api_key = os.getenv('OPENAI_API_KEY')  # set openai api key, stored via Pycharm settings "Build, Execution, Deployment"
cost_per_1000_input_tokens = 0.0015  # $0.0015 per 1000 input tokens, https://openai.com/pricing 2023-06-15 dollars
cost_per_1000_output_tokens = 0.002

# TODO:

# FIXME:

# ==========================
# Custom GPT Functions
# ==========================

### questions for the model and model function
def my_gpt_api_function(abstract_disease: str, mirna_oi: str):
    """
    Create a json schema
    https://community.openai.com/t/emulated-multi-function-calls-within-one-request/269582

    :param abstract_disease: The disease as a string
    :param mirna_oi: The microRNA of interest in the abstract
    """
    functions = [
        {
            "name": "mirco_rna_abstract_data",
            "description": "Retrieve study data from abstract",
            "parameters": {
                "type": "object",
                "properties": {
                    # Q1
                    "related_topic": {
                        "type": "string",
                        "enum": ["Yes", "No", "Not Sure"],
                        "description": f"Does the abstract deal with {abstract_disease} in humans?",
                    },
                    # Q2
                    "direction_upreg_downreg": {
                        "type": "string",
                        "enum": ["Upregulated", "Downregulated", "Not Given"],
                        "description": f"Is the given miRNA of interest {mirna_oi} up- (increased) or downregulated (decreased) in {abstract_disease}?",
                    },
                    # Q3
                    "primary_literature": {
                        "type": "string",
                        "enum": ["Yes", "No", "Not Sure"],
                        "description": "Is the abstract from an original research (no review, meta-analyses)?",
                    },
                    # Q4
                    "serum_plasma_tissue": {
                        "type": "string",
                        "enum": ["Serum", "Plasma", "Heart Tissue", "Not Sure"],
                        "description": f"Was human {mirna_oi} measured in serum, plasma, heart tissue or somewhere else?",
                    },
                    # Q5
                    "mortality": {
                        "type": "string",
                        "enum": ["Yes", "No", "Not Sure"],
                        "description": f"Is {mirna_oi} in the abstract clearly associated with mortality/death?",
                    },
                    # Q6
                    "measurement_type": {
                        "type": "string",
                        "enum": ["qPCR", "NGS", "Microarray", "Combination", "Not Sure"],
                        "description": f"Was human {mirna_oi} measured by qPCR, sequenced (NGS), used microarray?",
                    },
                    # Q7
                    "sample_size": {
                        "type": "integer",
                        "description": "What was the human sample size of the study?"
                    }
                },
                "required": ["related_topic", "direction_upreg_downreg", "primary_literature",
                             "serum_plasma_tissue", "mortality", "measurement_type", "sample_size"],
            },
        }
    ]
    return functions


# function to run model
# https://www.mlq.ai/gpt-function-calling-getting-started/

def run_conversation(GPT_MODEL, system_msg, user_msg, custom_gpt_function, retries=5):
    '''
    :param GPT_MODEL:           gpt model name
    :param system_msg:          describe the behavior of the AI assistant
    :param user_msg:            input text (abstract in our case)
    :param custom_gpt_function: gpt function calling syntax (structured questions here)
    :param retries:             number of retries  # could also use retry from tenacity library with the @retry decorator
    :return:                    response from gpt model
    '''

    try:  # could alternatively use requests library to call the API
        response = openai.ChatCompletion.create(
            model=GPT_MODEL,
            max_tokens=2048,
            temperature=0.7,
            top_p=1,
            messages=[{"role": "system", "content": system_msg},
                      {"role": "user", "content": user_msg}],
            functions=custom_gpt_function,
            function_call="auto",
        )
        # message = response["choices"][0]["message"]
        return response

    except Exception as e:
        if retries > 0:
            print(f"An error occurred: {str(e)}. {retries} left - Retrying in 20 seconds ...")
            time.sleep(20)
            print(f"Retrying {retries} more times")
            return run_conversation(GPT_MODEL, system_msg, user_msg, custom_gpt_function, retries - 1)
        else:  # after trying {retries} times
            print(f"An error occurred: {str(e)}. No retries left.", file=sys.stderr)
            return None


# ==========================
# Preparation
# ==========================

# %% load data
path2dat = "./data-literature/pubmed/"
list_of_files = sorted(os.listdir(path2dat))
concatenated_files = [path2dat + file for file in list_of_files if "-2023-06-16-" in file]

diseases = ["ACS", "CAD", "DCM", "HFrEF"]

# create pandas dataframe with disease abbreviations and full names
diseases_df = pd.DataFrame({"abbrev_disease": diseases,
                            "disease_name": ["Acute Coronary Syndrome", "Coronary Artery Disease",
                                             "Dilatative Cardiomyopathy", "Heart Failure"]
                            })

# ==========================
# Loop over diseases
# ==========================
# %%  loop over diseases and abstracts
current_date = datetime.now().strftime("%Y-%m-%d")

for i, disease in enumerate(diseases_df.disease_name.tolist()[1:], start=1):
    print(f"Processing abstracts for {disease} ({i+1}/{diseases_df.shape[0]})...")
    # ==========================
    df = pd.read_csv(concatenated_files[i], sep=";", low_memory=False)
    # initialize empty col names
    new_columns = ['usage_total_tokens', 'query_cost', 'time']
    # Using assign() method to add empty columns
    df = df.assign(**{col: pd.Series(dtype='float64') for col in new_columns})
    # Columns to change dtype to "object"
    columns_to_change = ['related_topic', 'direction_upreg_downreg', 'primary_literature', 'serum_plasma_tissue', 'mortality', 'measurement_type']
    df[columns_to_change] = df[columns_to_change].astype('object')
    df = df.drop('unrelated_topic01', axis=1)
    # df['time'] = df['time'].astype(float)
    # df.head()
    # df.info(); df.shape
    my_abstracts_grouped = df.groupby('PMID').agg({'miRNA': lambda x: list(x.unique())}).reset_index()
    # ==========================
    # loop over abstracts
    for index, row in df.iterrows():
        print(f"GPT is processing abstract and miRNA {index+1}/{df.shape[0]}...")
        # ==========================
        custom_gpt_function = my_gpt_api_function(abstract_disease=disease, mirna_oi=df["miRNA"][index])
        system_msg = 'You are a helpful assistant who understands data science and cardiovascular molecular biology.'
        #break
        # ==========================
        start_time = time.time()
        # run model
        gpt_response = run_conversation(
            GPT_MODEL=GPT_MODEL,
            system_msg=system_msg,
            user_msg=df.Abstract[index],
            custom_gpt_function=custom_gpt_function,
            retries=5
        )
        elapsed_time = time.time() - start_time  # in seconds
        # clean response
        message_dict = json.loads(gpt_response.choices[0]["message"]["function_call"]["arguments"])
        # test1 = pd.DataFrame.from_records([message_dict])  # convert dict to df
        prompt_tokens = gpt_response.usage["prompt_tokens"]
        completion_tokens = gpt_response.usage["completion_tokens"]
        total_tokens = prompt_tokens + completion_tokens
        query_cost = prompt_tokens/1000*cost_per_1000_input_tokens + completion_tokens/1000*cost_per_1000_output_tokens
        # ==========================
        # fill up dataframe
        # mirna data
        df.loc[index, 'related_topic'] = message_dict["related_topic"]
        df.loc[index, 'direction_upreg_downreg'] = message_dict["direction_upreg_downreg"]
        df.loc[index, 'primary_literature'] = message_dict["primary_literature"]
        df.loc[index, 'serum_plasma_tissue'] = message_dict["serum_plasma_tissue"]
        df.loc[index, 'mortality'] = message_dict["mortality"]
        df.loc[index, 'measurement_type'] = message_dict["measurement_type"]
        df.loc[index, 'sample_size'] = message_dict["sample_size"]
        # usage
        df.loc[index, 'usage_total_tokens'] = total_tokens
        df.loc[index, 'query_cost'] = query_cost
        df.loc[index, 'time'] = elapsed_time
        # ==========================
        # store original openai class object
        file_name_object = f"./output/gpt_response/{diseases[i]}/{current_date}-index{index}-gpt_response.pkl"
        with open(file_name_object, 'wb') as f:
            pickle.dump(gpt_response, f)
        # ==========================

    # store filled dataframe
    file_name_df_csv = f"./output/gpt_dataframe/{diseases[i]}/{current_date}-df-pubmed-gpt_response.csv"
    file_name_df_pkl = f"./output/gpt_dataframe/{diseases[i]}/{current_date}-df-pubmed-gpt_response.pkl"
    file_name_df_xlsx = f"./output/gpt_dataframe/{diseases[i]}/{current_date}-df-pubmed-gpt_response.xlsx"
    df.to_csv(file_name_df_csv, index=False)
    df.to_pickle(file_name_df_pkl)
    df.to_excel(file_name_df_xlsx, index=False)
    print(f"Finished processing abstracts for {disease}...\nThank you for your patience.\n==========================\n")
    # ========================== END OF LOOP ========================== #

