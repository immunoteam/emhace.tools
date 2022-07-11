#http://hla.alleles.org/nomenclature/naming.html
#HLA prefix, hyphen[-], gene[e.g. A, B, C, DRB1], separator[*], allele group (Field 1), field separator[:], specific HLA protein (Field 2)
#setwd("D:/CloudStation/netMHCpan-tools/")
#HLA-I
# allele = c("A0201","A020116","A0210102","A01:01","A01010101","A*01:01:01:01","A*02:01:01:02L","A*03:01:25","A*01:01:38L",
#            "A*02:210:01","A*02:01:100","A*03:01","A*03:10N","A*03:100","A*03:100N","A*11","HLA-A03:01","HLA-A1101","HLA-A*1101",
#            "HLA-A*110101","HLA-A13:01",
#            "DPA1*01:04","DPA1*01:21Q","DPA1*01:29N","DPA1*01:03:02","DPA1*01:03:01:01","DPA1*01:03:01:18Q","DPB1*08:01","DPB1*61:01N",
#            "DPB1*100:01","DPB1*120:01N","DPB1*697:01Q","DPB1*1009:01","DPB1*01:01:03","DPB1*1121:01N","DPB1*1148:01Q","DPB1*104:01:02",
#            "DPB1*786:01:01N","DPB1*01:01:01:01","DPB1*02:01:02:46Q","DPB1*04:01:01:24N","DPB1*104:01:01:01","HLA-DPA10103-DPB10101",
#            "HLA-DPA1*0103-DPB1*0101","HLA-DQA1*01:01/DQB1*05:01","HLA-DQA1*0101/DQB1*0501","HLA-DQA101:01/DQB105:01","HLA-DQA10101/DQB10501",
#            "DQA1*01:01/DQB1*05:01","HLA-DQA10103-DQB10101","HLA-DQA10103-DQB10101","HLA-DQA1*0103-DQB1*0101",
#            "DRB1*01:04", "DRB1*01:33N", "DRB1*01:91Q", "DRB1*01:100", "DRB1*03:156N", "DRB1*15:164Q", "DRB1*01:01:02", "DRB1*01:62:01N", 
#            "DRB1*03:100:01", "DRB1*01:01:01:01",
#            "HLA-G13:01","HLA-45:13FDA")


#HLA-II
#HLA-DPA10103-DPB10101 #netmhc

#DQ
#"DQA1*01:06", "DQA1*01:07Q", "DQA1*01:15N", "DQA1*01:01:03", "DQA1*01:01:01:01"
#"DQB1*05:04", "DQB1*05:100", "DQB1*05:110N", "DQB1*05:132Q", "DQB1*05:01:02", "DQB1*06:118:01", "DQB1*05:01:01:01", "DQB1*03:01:01:21Q", "DQB1*03:263:01:01"
#HLA-DQA10103-DQB10101 #netmhc

#DR
#DRB1_0101 #netmhc



#Output format types
#2: HLA-A*01:01 - 4digit
#3: HLA-A01:01 - netmhc
#4: A0101 - simple

# packs = c("fastmatch")
# invisible(lapply(packs, require, character.only = TRUE))

###########
#"http://www.cbs.dtu.dk/services/NetMHCIIpan/alleles_name.list"
#"http://www.cbs.dtu.dk/services/NetMHCpan/MHC_allele_names.txt"

SimplifyAlleles = function(allele, output.format = "netmhc") {
  
  packs <- c("fastmatch","stringr")
  invisible(lapply(packs, require, character.only = TRUE))
  netmhc_allowed_alleles = readLines("/home/workstation/netMHCpan-tools/netmhcpan_allowed_alleles.txt")
  mysplit = function(s) {unlist(strsplit(s,":"),use.names = F)}
  
  sapply(allele, function(a) {
    tempallele = gsub("HLA-|\\*|_|L|N","",a)
    if(substr(tempallele,1,1) == "D") {
      hlaclass = "ii"
    } else if(substr(tempallele,1,1) %in% c("A","B","C")) {
      hlaclass = "i"
    } else {
      hlaclass = "unknown"
    }
    
    if(hlaclass == "unknown") {
      tempallele = "invalid or incomplete allele"
    } else if(hlaclass == "i") {
      #HLA-I
      tempallele = gsub("/|-","",tempallele)
      if(str_count(string = tempallele, pattern = "[ABC]") != 1 | length(intersect(unlist(strsplit(tempallele,""),use.names = F), setdiff(LETTERS,c("A","B","C")))) > 0) {
        tempallele = "invalid or incomplete allele"
      } else if(output.format == "simple") {
        if(grepl(":", tempallele)) tempallele = paste0(mysplit(tempallele),collapse = "")
      } else if(output.format == "4digit") {
        if(grepl(":", tempallele)) {
          tempallele = paste0(substr(tempallele,1,3),mysplit(tempallele)[2])
        } else if(nchar(tempallele)%%2 == 0) {
          tempallele = substr(tempallele,1,6)
        } else {
          tempallele = substr(tempallele,1,5)
        }
      } else if(output.format == "netmhc") {
        if(grepl(":", tempallele)) {
          tempallele = paste0("HLA-", paste0(mysplit(tempallele)[1:2],collapse = ":"))
        } else if(nchar(tempallele)%%2 == 0) {
          tempallele = paste0("HLA-", substr(tempallele,1,3), ":", substr(tempallele,4,6))
        } else {
          tempallele = paste0("HLA-", substr(tempallele,1,3), ":", substr(tempallele,4,5))
        }
        if(is.na(fmatch(tempallele,netmhc_allowed_alleles))) tempallele = "NetMHCpan is unable to predict binding to this allele."
      }

    } else {
      #HLA-II
      
      
      if(length(intersect(unlist(strsplit(tempallele,""),use.names = F), setdiff(LETTERS,c("A","B","D","P","Q","R")))) > 0) {
        tempallele = "invalid allele"
      } else if(substr(tempallele,1,3) == "DRB") {
        tempallele = gsub("/|-","",tempallele)
        
        intersect(unlist(strsplit(tempallele,""),use.names = F), setdiff(LETTERS,c("A","B","D","P","Q","R")))
        
        if(output.format == "simple") {
          if(grepl(":", tempallele)) tempallele = paste0(mysplit(tempallele),collapse = "")
        } else if(output.format == "4digit") {
          if(grepl(":", tempallele)) {
            tempallele = paste0(substr(tempallele,1,6),mysplit(tempallele)[2])
          } else if(nchar(tempallele)%%2 == 0) {
            tempallele = substr(tempallele,1,8)
          } else {
            tempallele = substr(tempallele,1,9)
          }
        }  else if(output.format == "netmhc") {
          if(grepl(":", tempallele)) {
            tempallele = paste0("DRB1_",substr(tempallele,5,6),mysplit(tempallele)[2])
          } else if(nchar(tempallele)%%2 == 0) {
            tempallele = paste0("DRB1_",substr(tempallele,5,8))
          } else {
            tempallele = paste0("DRB1_",substr(tempallele,5,9))
          }
        }
        if(str_count(tempallele, "D") != 1 | str_count(string = tempallele, pattern = "R") != 1 | str_count(string = tempallele, pattern = "B") != 1 | length(intersect(unlist(strsplit(tempallele,""),use.names = F), setdiff(LETTERS,c("A","B","D","P","Q","R")))) > 0 | str_count(tempallele, "[A-Z]")>3) {
          tempallele = "invalid or incomplete allele"
        }
      } else if((grepl("DPA", tempallele) & grepl("DPB", tempallele)) | (grepl("DQA", tempallele) & grepl("DQB", tempallele))) {
        output.format = "netmhc"
        if(grepl("/",tempallele)) {
          tempallele1 = strsplit(tempallele,"/")[[1]][1]
          tempallele2 = strsplit(tempallele,"/")[[1]][2]
        } else if(grepl("-",tempallele)) {
          tempallele1 = strsplit(tempallele,"-")[[1]][1]
          tempallele2 = strsplit(tempallele,"-")[[1]][2]
        } else {
          tempallele1 = paste0("D",strsplit(tempallele,"D")[[1]][2])
          tempallele2 = paste0("D",strsplit(tempallele,"D")[[1]][3])
        }
        #netmhc_check
        if(grepl(":", tempallele1)) {
          tempallele1_netmhc = paste0(substr(tempallele1,1,4), "*", substr(tempallele1,5,nchar(tempallele1)))
        } else {
          tempallele1_netmhc = paste0(substr(tempallele1,1,4), "*", substr(tempallele1,5,6), ":", substr(tempallele1,7,nchar(tempallele1)))
        }
        if(grepl(":", tempallele2)) {
          tempallele2_netmhc = paste0(substr(tempallele2,1,4), "*", substr(tempallele2,5,nchar(tempallele1)))
        } else {
          tempallele2_netmhc = paste0(substr(tempallele2,1,4), "*", substr(tempallele2,5,6), ":", substr(tempallele2,7,nchar(tempallele1)))
        }
        if(grepl(":", tempallele1)) tempallele1 = paste0(mysplit(tempallele1),collapse = "")
        if(grepl(":", tempallele2)) tempallele2 = paste0(mysplit(tempallele2),collapse = "")
        tempallele = paste("HLA", tempallele1, tempallele2, sep = "-")
        if(is.na(fmatch(tempallele1_netmhc,netmhc_allowed_alleles)) | is.na(fmatch(tempallele2_netmhc,netmhc_allowed_alleles))) tempallele = "NetMHCpan is unable to predict binding to this allele OR invalid allele."
      } else {
        tempallele = "invalid or incomplete allele"
      }
    }
    #if(nchar(tempallele)<5) print("Allele list contains allele with low resolution (2-digit)")
    if((output.format %in% c("simple","4digit","netmhc")) == FALSE) tempallele = "invalid output format"
    #if(grepl("DQ", tempallele) | grepl("DP", tempallele)) cat("output format: netmhc")
    return(tempallele)
  })
}

# SimplifyAlleles(allele, output.format = "netmhc")
# 
# #test
# allele = c("sldkajsdlk0", "A0201:054512", "A0201:054512A", "A0201:054512D", "AB0201:054512", "A0201:054512HLA", "DQSA0121")
#            
# 
# library(XML)
# hlaii = readHTMLTable("http://hla.alleles.org/alleles/class2.html")
# View(hlaii$table1)
# table(nchar(hlaii$table1$`HLA-DRB1`))
# View(hlaii$table1$`HLA-DRB1`[nchar(hlaii$table1$`HLA-DRB1`) == 10])
# table(nchar(hlaii$table1$`HLA-DQA1`))
# 
