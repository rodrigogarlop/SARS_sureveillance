# First, load the complete table
library("readxl")
# library("xtable")
xlsx <- "01_data/2023_03_16_DB_MexCov2.xlsx"
df <- as.data.frame(read_excel(xlsx, sheet = 1))
# df <- read.table("01_data/2022_04_27_DB_MexCov2.tsv",header=T, sep='\t', skip=0, comment.char='',fill=T, check.names=FALSE, stringsAsFactors = FALSE)
# df <- read.table(unz("01_data/17Enero2022_Metadata_BDMexCov2_FINAL.zip", "17Enero2022_Metadata_BDMexCov2_FINAL.tsv"), header=T, quote="\"", sep="\t", check.names=FALSE, fill=TRUE)
print("Starting dimensions:")
dim(df)
# [1] 86468    33
# forDel <- sort(c(grep("aboratorio", colnames(df)), grep("Mutaciones",colnames(df)), grep("Deleciones",colnames(df)),grep("Inserciones",colnames(df)), grep("Editar", colnames(df)))) # If we want to remove mutations (to reduce table size)
forDel <- sort(c(grep("aboratorio", colnames(df)), grep("Editar", colnames(df)))) # Alt ver, do not remove mutations
df <- df[,-forDel] # remove info on mutations (as some of these have faulty characters that break some downstream analyses)
# max <- grep("Delta_Sub",colnames(df))-1 # All columns starting from Delta_Sub should be removed as they are older versions and will be recalculated in subsequent steps
# df <- df[,1:max]
dim(df)
# [1] 86468    28
# Remove rare punctuations or accents
strip <- function(string){ # Takes a string object and strips rare characters
	string <- iconv(string,from="UTF-8",to="ASCII//TRANSLIT") # Remove accents and rare chars
# 	string <- tolower(string) # change to lower
	string <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", string, perl=TRUE) # Remove multiple spaces
	return(string)
}
names(df) <- strip(names(df))
df[,"Estado"] <- strip(df[,"Estado"])
df[,"Municipio"] <- strip(df[,"Municipio"])
df[,"Edad"] <- strip(df[,"Edad"])
# Fix the Ambulatory column
df[df[,"Estatus del paciente"]=="Asymptomatic","Estatus del paciente"] <- "Ambulatory"
df[df[,"Estatus del paciente"]=="Asymptomatic and Ambulatory","Estatus del paciente"] <- "Ambulatory"
df[df[,"Estatus del paciente"]=="Outpatient","Estatus del paciente"] <- "Ambulatory"
df[df[,"Estatus del paciente"]=="Mild","Estatus del paciente"] <- "Ambulatory"
df[df[,"Estatus del paciente"]=="Symptomatic and Ambulatory","Estatus del paciente"] <- "Ambulatory"
df[df[,"Estatus del paciente"]=="Symptomatic and Hospitalized","Estatus del paciente"] <- "Hospitalized"
df[df[,"Estatus del paciente"]=="Released","Estatus del paciente"] <- "Hospitalized"
df[df[,"Estatus del paciente"]=="Fatal","Estatus del paciente"] <- "Deceased"
df[df[,"Estatus del paciente"]=="Live","Estatus del paciente"] <- "unknown"
df[df[,"Estatus del paciente"]=="Moderate","Estatus del paciente"] <- "unknown"
df[df[,"Estatus del paciente"]=="Symptomatic","Estatus del paciente"] <- "unknown"
# # # # Add a new column to emulate the old input table for recycling the script Filter_Delta.R
# # # df["FLAG"] <- df[,"Linaje Pangolin"]=="B.1.617.2"
# # # df[grep("^AY\\.",df[,"Linaje Pangolin"]),"FLAG"] <- TRUE
# # # df[df[,"FLAG"],"ID"] <-  df[df[,"FLAG"],"Accession ID"]
# # # df[df[,"FLAG"],"New_lineage"] <-  df[df[,"FLAG"],"Linaje Pangolin"]
# # # df[is.na(df)] <- ""
# # # df <- df[,-which(names(df)=="FLAG")]
print("Ending dimensions:")
dim(df)
# [1] 86468    28
write.table(df, "01_data/preFiltered.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
