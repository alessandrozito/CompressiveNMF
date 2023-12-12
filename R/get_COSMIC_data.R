# This file scrapes the data from the COSMIC database. It imports bothe the signatures and 
# the details on the proposed aetiology to faciltate interpretation.
library(rvest) # to scrape the data
library(tidyverse)

################################################################################
# Source the data - COSMIC v3.4
################################################################################ 

# Dowload the cosmic database
#'' @example cosmic_data <- import_cosmic_data()
import_cosmic_data <- function(url = "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.4_SBS_GRCh37.txt") {
  # Read the tab-delimited text file from the provided URL
  cosmic_data <- as.data.frame(read.table(url, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
  # Return the imported data
  channel <- cosmic_data$Type
  channel_short <- str_replace_all(str_extract(channel, "\\[(.*?)\\]"), "\\[|\\]", "")
  channel_new <- str_replace(channel, "\\[.*?\\]", function(match) {
    first_letter_after_open_bracket <- str_sub(match, start = 2, end = 2)
    return(first_letter_after_open_bracket)
  })
  cosmic_data <-  cbind("Channel" = channel, "Mutation" = channel_short, "Triplet" =channel_new , cosmic_data[, -1])
  return(cosmic_data)
}


# There are 87 signatures. What is important now is to filter out the ones that are not necessary
# Take the relevant one like SBS5. 

# To understand how it works, we need to scrape some data from the internet and 
# build a table
scrape_Cosmic <- function(){
  cosmic_data <- import_cosmic_data()
  singature_names <- colnames(cosmic_data)[-1]
  url <- "https://cancer.sanger.ac.uk/signatures/sbs/"
  
  data_sign <- data.frame()
  for(i in 1:length(singature_names)){
    # get signature name
    sign <- singature_names[i]
    url_sign <- paste0(url, str_to_lower(sign), "/")
    # Scrape the comments
    webpage <- read_html(url_sign)
    # Get etiology
    aet <- webpage %>%
      html_elements("#proposed-aetiology") %>% 
      html_elements("p") %>% 
      html_text2() %>%
      matrix(ncol = 2, byrow = TRUE, dimnames = list(NULL,c("Proposed_aetiology", "Comments"))) %>%
      as.data.frame()
    # Get the table info 
    tab <- webpage %>% 
      html_node("table") %>% 
      html_table() %>%
      as.data.frame()
    tab_full <- data.frame("signature" = sign,
                           "aetiology_short" = c(tab[7,2]),
                           "id_study" = tab[1,2], 
                           "COSMIC_version" = tab[1,4],
                           "Validated" = tab[5,2],
                           "Replicated" = tab[5,3],
                           "Experimental_validation" = tab[9,2])
    data_sign <- rbind(data_sign, cbind(tab_full, aet)) 
  }
  return(data_sign)
}

#"#022B7B"

#aetiology_data <- data_aet
#save(x = aetiology_data, file = "data/COSMICsignatures_aetiology.rdata")





