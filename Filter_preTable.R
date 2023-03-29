# Started 2023-02-03 by Rodrigo Garcia-Lopez for Prof. Arias laboratory at IBt, UNAM
# Under GNU GPLv3 license
# This script takes a GISAID-based metadata table that has lineages and appends a table with
# LOAD FUNCTIONS
get_full_lineage <- function(s_org,dict){
	for(i in 1:nrow(dict)){
# 		print(cat(s_org, paste0("^",dict[i,1],"\\."),dict[i,2],"\n"))
		s_org <- sub("Unassigned", "unassigned", s_org)
		s_org <- sub(paste0("^",dict[i,1]), dict[i,2], s_org)
	}
	return(s_org)
}
# LOAD INPUTS
df <- read.table("01_data/preFiltered.tsv", header=T, quote="", sep="\t", check.names=FALSE)
# df <- read.table("01_data/Xai_tabla_filtrada_Delta_fixed5.tsv",header=T, sep='\t', skip=0, comment.char='',fill=FALSE, check.names=FALSE, stringsAsFactors = FALSE)
full_lineages <- read.table("01_data/variant_full_lineage.tsv",header=FALSE, sep='\t', skip=0, comment.char="",fill=T, check.names=FALSE, stringsAsFactors = FALSE, row.names=NULL) # Load dictionary (c1=alias, c2=full_lineage, no headers)
print("Starting total items:")
dim(df)
# [1] 86468    28
dim(full_lineages) # Items in the dictionary
# [1] 128   2
targetCol <- grep("lin", tolower(colnames(df))) # Get the target column (change string accordingly)
df[,"New_Lineage"] <- get_full_lineage(df[,targetCol], full_lineages)
# UPDATE 2023-03-22: Now, append additional information for sequences flagged as "Under investigation". As of this date, there are 190 such items but only 188 in our data. There may be marked due to different reasons including having multiple lineages (potential recombinants, marked ML), having frame shifts (FS), having submission dates earlier than the actual sample collection dates (CD), and truncated spike proteins (ST).
UnderInv <- read.table("01_data/2023-03-22_Mex_UnderInvestigation_with_reason_e.g.ML.tsv", header=T, quote="", sep="\t", check.names=FALSE)
UnderInv <- UnderInv[,c("Accession ID" ,"Under investigation")] # Only keep the useful columns (epi id and Under Investigation cols)
rownames(UnderInv) <- UnderInv[,1] # Use as dictionary
UnderInv <- UnderInv[UnderInv[,2]=="ML",] # KEEP ONLY MLs
df[,"ML"] <- UnderInv[df[,"Accession ID"],2] # Append a column for those with actual info
temp <- df[,"ML"]=="ML";temp[is.na(temp)] <- FALSE # Get a vector with ML positions
df[temp,"Linaje Pangolin"] <- "MultLin" # Now, mark those with Multiple Lineages
df[temp,"New_Lineage"] <- "MultLin" # Now, mark those with Multiple Lineages

# df[!is.na(df[,"ML"]),"ML_flag"] <- "ML" # Add a flag
write.table(df, "01_data/preFiltered_full_lineage.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

df[grepl("^B\\.1\\.617\\.2$",df[,"New_Lineage"]),"Delta"] <- "Delta" # Flag those that are delta variants
df[grepl("^B\\.1\\.617\\.2\\.",df[,"New_Lineage"]),"Delta"] <- "Delta"
df[is.na(df[,"Delta"]),"Delta"] <- "Others"
# table(df[,"Delta"])
#   Delta Others
#  25706  60762
# Create subcategories
df[,"Delta_sub"] <- df[,"Delta"]
# sort(table(df[df[,"Delta_sub"]=="Delta","Linaje Pangolin"]),decreasing=T)
#      AY.20      AY.26     AY.100       AY.3     AY.113     AY.103      AY.44
#      11620       5605       3235       1543        980        595        315
#      AY.62  B.1.617.2     AY.122      AY.25    AY.25.1      AY.39      AY.43
#        270        246        216        163        161        160        110
#      AY.13      AY.47     AY.118      AY.75   AY.119.2     AY.119       AY.4
#         84         52         33         25         21         19         17
#      AY.15   AY.122.4      AY.42      AY.74     AY.127      AY.38     AY.117
#         16         15         15         13         12         12         11
#      AY.14      AY.54   AY.119.1       AY.2     AY.114     AY.121    AY.39.2
#         11         11          9          9          8          8          7
#     AY.3.1      AY.53       AY.5    AY.98.1   AY.116.1     AY.101     AY.126
#          6          6          5          5          4          3          3
#      AY.33    AY.99.2     AY.105   AY.122.6      AY.34      AY.35     AY.4.2
#          3          3          2          2          2          2          2
#     AY.4.9      AY.46    AY.46.4      AY.48      AY.52   AY.120.1   AY.121.1
#          2          2          2          2          2          1          1
#     AY.124 AY.124.1.1     AY.125      AY.27     AY.3.3      AY.32    AY.34.1
#          1          1          1          1          1          1          1
#    AY.34.2      AY.37    AY.39.1    AY.4.10       AY.6      AY.64      AY.66
#          1          1          1          1          1          1          1
#      AY.67      AY.73      AY.77     AY.9.2      AY.98
#          1          1          1          1          1
df[df[,"Delta_sub"]=="Delta","Delta_sub"] <- "Other Delta" # First, change the name
df[df[,"Linaje Pangolin"]=="AY.20","Delta_sub"] <- "AY.20"
df[df[,"Linaje Pangolin"]=="AY.26","Delta_sub"] <- "AY.26"
df[df[,"Linaje Pangolin"]=="AY.3","Delta_sub"] <- "AY.3"
df[df[,"Linaje Pangolin"]=="AY.113","Delta_sub"] <- "AY.113"
df[df[,"Linaje Pangolin"]=="AY.103","Delta_sub"] <- "AY.103"
df[df[,"Linaje Pangolin"]=="AY.100","Delta_sub"] <- "AY.100"
df[df[,"Linaje Pangolin"]=="AY.4","Delta_sub"] <- "AY.44"
df[df[,"Linaje Pangolin"]=="AY.62","Delta_sub"] <- "AY.62"
df[df[,"Linaje Pangolin"]=="B.1.617.2","Delta_sub"] <- "B.1.617.2"
df[df[,"Linaje Pangolin"]=="AY.122","Delta_sub"] <- "AY.122"
df[df[,"Linaje Pangolin"]=="AY.39","Delta_sub"] <- "AY.39"
df[df[,"Linaje Pangolin"]=="AY.25","Delta_sub"] <- "AY.25"
df[df[,"Linaje Pangolin"]=="AY.25.1","Delta_sub"] <- "AY.25" # This is joined with AY.25
df[df[,"Linaje Pangolin"]=="AY.43","Delta_sub"] <- "AY.43"

# sort(table(df[,"Delta_sub"]), decreasing=TRUE)
#      Others       AY.20       AY.26      AY.100        AY.3      AY.113
#       60762       11620        5605        3235        1543         980
# Other Delta      AY.103       AY.25       AY.62   B.1.617.2      AY.122
#         785         595         324         270         246         216
#       AY.39       AY.43       AY.44
#         160         110          17

# Create a column for tracking omicron variants
df[grepl("^B\\.1\\.1\\.529$",df[,"New_Lineage"]),"Omicron"] <- "Omicron" # Flag those that are omicron variants
df[grepl("^B\\.1\\.1\\.529\\.",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
# table(df[grepl("^X",df[,"New_Lineage"]),"Linaje Pangolin"]) # Recombinants are checked manually, as most are Omicron
#       XAF       XAG       XAH       XAJ       XAM       XAP       XAS       XAZ
#         2         2         4         1        12         4        78         2
#        XB       XBB     XBB.1   XBB.1.1   XBB.1.2   XBB.1.4 XBB.1.4.1   XBB.1.5
#       326        15       607         4        58        22         1       595
#   XBB.1.6   XBB.1.9     XBB.2   XBB.2.2     XBB.3       XBD       XBF        XZ
#        12        11        43         4         4         1         2         1

df[grepl("^X",df[,"New_Lineage"]),"Recomb"] <- "Recomb" # add a new column for recombinants (all which start with X)
df[is.na(df[,"Recomb"]),"Recomb"] <- "Others" # the rest are tagged as "Others"
# and append the ML if present
df[temp,"Recomb"] <- "MultLin"
# table(df[,"Recomb"])
#     ML Others Recomb
#     89  84572   1807
df[,"Recomb_two"] <- df[,"Recomb"] # add a new column to get details on recombinants
df[grepl("^X",df[,"New_Lineage"]),"Recomb_two"] <- df[grepl("^X",df[,"New_Lineage"]),"New_Lineage"]
df[temp,"Recomb_two"] <- "MultLin"

# table(df[,"Recomb_two"])
#        ML    Others       XAF       XAG       XAH       XAJ       XAM       XAP
#        89     84572         2         1         4         1         9         4
#       XAS       XAZ        XB       XBB     XBB.1   XBB.1.1   XBB.1.2   XBB.1.4
#        78         2       326        15       607         4        58        22
# XBB.1.4.1   XBB.1.5   XBB.1.6   XBB.1.9     XBB.2   XBB.2.2     XBB.3       XBD
#         1       595        12        11        43         4         4         1
#       XBF        XZ
#         2         1

# Now, add known Omicron recombinants to the omicron List (as of 2023-02-07, no delta-omicron recombintants have been reported
df[grepl("^XZ$",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
df[grepl("^XAF$",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
df[grepl("^XAG$",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
df[grepl("^XAH$",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
df[grepl("^XAJ$",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
df[grepl("^XAM$",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
df[grepl("^XAP$",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
df[grepl("^XAS$",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
df[grepl("^XAZ$",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
df[grepl("^XBB",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
df[grepl("^XBD$",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
df[grepl("^XBF$",df[,"New_Lineage"]),"Omicron"] <- "Omicron"
df[is.na(df[,"Omicron"]),"Omicron"] <- "Others"
# table(df[,"Omicron"])
# Omicron  Others
#   39960   46508

# sort(table(df[df[,"Omicron"]=="Omicron","Linaje Pangolin"]),decreasing=T)
#    BA.1.1 BA.2.12.1    BA.5.1   BA.1.15      BA.1  BA.5.1.6      BA.2  BA.5.2.1
#      9007      3052      2669      2505      2054      1492      1424      1412
#    BA.2.9  BA.5.2.9    BQ.1.1 BA.5.1.23     XBB.1   XBB.1.5  BA.5.1.1    BA.4.1
#      1237      1133       964       650       607       595       574       519
#    BW.1.1    BA.5.2      BW.1   BQ.1.14    BA.4.2      BQ.1 BA.5.1.25    BA.4.4
#       478       455       426       386       350       340       335       287
#      BF.1 BA.5.1.30    BA.5.6   BA.1.17      BU.1  BA.4.1.8    BA.5.5 BA.1.1.18
#       279       232       215       199       194       192       192       190
#   BA.1.20    BA.2.3    BE.1.4     BF.10      BE.1      BF.8 BQ.1.1.10 BQ.1.1.13
#       189       180       177       154       148       147       146       144
# BA.5.1.22  BA.5.1.3      BE.2 BA.1.17.2 BA.5.2.22      BA.4      BE.3    BQ.1.2
#       132       120       120       117       116       112       106       106
#   BA.2.18    BA.4.6       XAS  BQ.1.1.4 BA.1.15.2   BA.1.18      BA.5   BA.1.21
#        94        87        78        72        68        66        65        63
# BQ.1.1.18   BQ.1.10     BF.26   XBB.1.2    BE.1.1   BQ.1.12  BA.1.1.1     XBB.2
#        63        61        59        58        53        51        43        43
#      BF.7  BQ.1.1.5 B.1.1.529      CK.1     BF.27   BA.1.13      BN.1 BA.2.10.1
#        41        40        39        39        38        37        36        35
#      BF.5    BQ.1.3  BA.5.6.2    BA.5.3  BN.1.3.1 BA.5.2.23 BA.1.1.14 BA.5.1.10
#        34        34        33        32        32        30        29        28
# BA.5.2.21    BN.1.3   BA.2.65     BF.14   BQ.1.13    CM.8.1  BQ.1.1.3    BQ.1.5
#        28        28        27        27        27        27        26        26
#      DE.1   BA.2.13  BA.5.1.4 BA.1.1.16 BA.5.1.15     BF.21    CH.1.1  BA.5.3.1
#        26        25        25        24        24        24        24        22
#    BQ.1.8  CK.2.1.1   XBB.1.4 BA.5.2.20 BA.2.40.1  BA.4.6.5  BA.2.9.3   BQ.1.25
#        22        22        22        21        20        20        19        19
#    BN.3.1   BA.2.21    BN.1.4     BF.28 BQ.1.1.15   BA.2.10 BQ.1.1.32       XBB
#        18        17        17        16        16        15        15        15
#  BA.1.1.2   BA.2.56  BA.4.1.1 BA.5.1.12  BE.1.4.3  BQ.1.1.1   BQ.1.11   BA.2.22
#        14        14        14        14        14        14        14        13
#      BG.5    BA.1.5  BA.5.2.6      BE.5       XAM   XBB.1.6  BA.2.3.6   BA.2.36
#        13        12        12        12        12        12        11        11
# BA.5.1.21   XBB.1.9    BA.2.1 BA.5.1.24     BE.10    BN.1.5    BQ.1.6   BA.2.12
#        11        11        10        10        10        10        10         9
#  BA.4.1.6  BF.7.4.1      DB.2  BA.5.1.5  BA.5.2.3    BA.5.8    BN.1.2      CA.7
#         9         9         9         8         8         8         8         8
#      CN.1 BA.2.13.1   BA.2.23   BA.2.38   BA.2.71 BA.5.2.34  BE.1.1.2    BF.7.5
#         8         7         7         7         7         7         7         7
#   BQ.1.15 BA.2.3.10 BA.2.75.2 BA.5.1.19 BA.5.2.31  BA.5.2.8  BE.1.1.1  BE.1.2.1
#         7         6         6         6         6         6         6         6
#      BF.4 BQ.1.1.27 BQ.1.13.1  BA.5.1.2   BA.5.10     BF.31  BQ.1.1.2  BQ.1.1.7
#         6         6         6         5         5         5         5         5
# BA.1.1.10 BA.1.15.1 BA.1.15.3    BA.2.6   BA.2.66 BA.5.2.13 BA.5.2.24 BA.5.2.33
#         4         4         4         4         4         4         4         4
#  BA.5.5.1  BA.5.5.3   BF.10.1    BF.7.4    BN.1.9      CT.1       XAH       XAP
#         4         4         4         4         4         4         4         4
#   XBB.1.1   XBB.2.2     XBB.3   BA.1.19 BA.2.3.15  BA.2.3.4   BA.2.37   BA.2.72
#         4         4         4         3         3         3         3         3
#  BA.2.9.7 BA.5.1.27 BA.5.1.28 BA.5.2.28  BA.5.5.2      BE.4      BE.9     BF.25
#         3         3         3         3         3         3         3         3
#      BG.2      BG.4      BK.1    BM.1.1 BQ.1.1.26   BQ.1.18   BQ.1.19   BQ.1.23
#         3         3         3         3         3         3         3         3
#      BR.2      BY.1  BY.1.1.1      CF.1      DF.1      DM.1 BA.1.1.15  BA.1.1.6
#         3         3         3         3         3         3         2         2
#   BA.2.26   BA.2.31   BA.2.44   BA.2.52 BA.5.10.1 BA.5.2.16 BA.5.2.35 BA.5.2.36
#         2         2         2         2         2         2         2         2
#  BA.5.6.1    BA.5.9    BE.1.2    BF.1.1   BF.11.2     BF.12     BF.22      BG.1
#         2         2         2         2         2         2         2         2
#      BG.3 BQ.1.1.11 BQ.1.1.23 BQ.1.1.28  BQ.1.1.6  BQ.1.1.8 BQ.1.10.1      CA.1
#         2         2         2         2         2         2         2         2
#  CH.1.1.2      CK.3      CL.1    CM.5.2    CR.1.1      CV.2    DF.1.1      DL.1
#         2         2         2         2         2         2         2         2
#       XAF       XAG       XAZ       XBF  BA.1.1.7  BA.1.1.9   BA.1.10   BA.1.14
#         2         2         2         2         1         1         1         1
# BA.1.14.1 BA.1.14.2   BA.1.16    BA.1.6    BA.1.7    BA.2.2 BA.2.23.1   BA.2.24
#         1         1         1         1         1         1         1         1
# BA.2.3.11 BA.2.3.17  BA.2.3.2   BA.2.35    BA.2.4   BA.2.49    BA.2.5   BA.2.50
#         1         1         1         1         1         1         1         1
#   BA.2.51   BA.2.57    BA.2.7   BA.2.75 BA.2.75.6   BA.2.78    BA.2.8   BA.2.81
#         1         1         1         1         1         1         1         1
#   BA.2.82  BA.2.9.5  BA.2.9.6  BA.4.6.1    BA.4.7 BA.5.1.17 BA.5.1.18 BA.5.1.20
#         1         1         1         1         1         1         1         1
#  BA.5.1.8   BA.5.11 BA.5.2.26 BA.5.2.29 BA.5.2.32 BA.5.2.37 BA.5.2.39 BA.5.2.47
#         1         1         1         1         1         1         1         1
#  BA.5.6.3    BE.1.3  BE.1.4.1  BE.1.4.2     BF.11   BF.11.1     BF.13      BF.2
#         1         1         1         1         1         1         1         1
#     BF.32    BF.7.3  BF.7.4.2    BF.7.6    BF.7.7      BL.6      BM.1  BM.4.1.1
#         1         1         1         1         1         1         1         1
#    BN.1.7 BQ.1.1.29 BQ.1.1.31      BR.1      BS.1    CA.3.1      CD.2  CH.1.1.1
#         1         1         1         1         1         1         1         1
#    CK.2.1    CM.4.1      CQ.2      CR.1      CV.1      DB.1      DG.1      DJ.1
#         1         1         1         1         1         1         1         1
#    DJ.1.1       XAJ XBB.1.4.1       XBD        XZ
#         1         1         1         1         1

df[,"Omicron_sub"] <- df[,"Omicron"]
df[df[,"Omicron_sub"]=="Omicron","Omicron_sub"] <- df[df[,"Omicron_sub"]=="Omicron","Linaje Pangolin"]
df[,"Omicron_two"] <- df[,"Omicron"]
df[grep("B\\.1\\.1\\.529$",df[,"New_Lineage"]), "Omicron_two"] <- "BA.1.x"
df[grep("B\\.1\\.1\\.529\\.1$",df[,"New_Lineage"]), "Omicron_two"] <- "BA.1.x"
df[grep("B\\.1\\.1\\.529\\.1\\.",df[,"New_Lineage"]), "Omicron_two"] <- "BA.1.x"
df[grep("B\\.1\\.1\\.529\\.2$",df[,"New_Lineage"]), "Omicron_two"] <- "BA.2.x"
df[grep("B\\.1\\.1\\.529\\.2\\.",df[,"New_Lineage"]), "Omicron_two"] <- "BA.2.x"
df[grep("B\\.1\\.1\\.529\\.3$",df[,"New_Lineage"]), "Omicron_two"] <- "BA.3.x"
df[grep("B\\.1\\.1\\.529\\.3\\.",df[,"New_Lineage"]), "Omicron_two"] <- "BA.3.x"
df[grep("B\\.1\\.1\\.529\\.4$",df[,"New_Lineage"]), "Omicron_two"] <- "BA.4.x"
df[grep("B\\.1\\.1\\.529\\.4\\.",df[,"New_Lineage"]), "Omicron_two"] <- "BA.4.x"
df[grep("B\\.1\\.1\\.529\\.5$",df[,"New_Lineage"]), "Omicron_two"] <- "BA.5.x"
df[grep("B\\.1\\.1\\.529\\.5\\.",df[,"New_Lineage"]), "Omicron_two"] <- "BA.5.x"
df[as.logical((df[,"Omicron_two"]=="Omicron")*(df[,"Recomb"]=="Recomb")),"Omicron_two"] <- "Omicron_Recomb"
# sort(table(df[,"Omicron_two"]))
# Omicron_Recomb         BA.4.x         BA.2.x         BA.1.x         BA.5.x
#           1485           1592           6571          14684          15628
#         Others
#          46508

df[,"Omicron_three"] <- df[,"Omicron_two"] # This will get some of the most abundant and split BA.5 specifically
df[df[,"Omicron_three"]=="BA.3.x","Omicron_three"] <- "Other Omicron"
df[grep("B\\.1\\.1\\.529$",df[,"New_Lineage"]), "Omicron_three"] <- "Other Omicron"
df[grep("B\\.1\\.1\\.529\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.1"
df[grep("B\\.1\\.1\\.529\\.1\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.1.1"
# df[grep("B\\.1\\.1\\.529\\.1\\.1\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BA.1.1.x"
df[grep("B\\.1\\.1\\.529\\.1\\.15$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.1.15"
# df[grep("B\\.1\\.1\\.529\\.1\\.15\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BA.1.15.x"
df[grep("B\\.1\\.1\\.529\\.2$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.2"
df[grep("B\\.1\\.1\\.529\\.2\\.9$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.2.9"
# df[grep("B\\.1\\.1\\.529\\.2\\.9\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BA.2.9.x"
df[grep("B\\.1\\.1\\.529\\.2\\.12\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.2.12.1"
# df[grep("B\\.1\\.1\\.529\\.2\\.12\\.1\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BA.2.12.1.x"
# df[grep("B\\.1\\.1\\.529\\.2\\.75$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.2.75"
# df[grep("B\\.1\\.1\\.529\\.2\\.75\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BA.2.75.x"
df[grep("B\\.1\\.1\\.529\\.4$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.4/BA.4.x"
# df[grep("B\\.1\\.1\\.529\\.4\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.4.1"
# df[grep("B\\.1\\.1\\.529\\.4\\.1\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BA.4.1.x" # This was the only BA.4 lineage with several cases
# df[grep("B\\.1\\.1\\.529\\.4\\.2$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.4.2"
# df[grep("B\\.1\\.1\\.529\\.4\\.4$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.4.4"
df[grep("B\\.1\\.1\\.529\\.5$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5/BA.5.x"
df[grep("B\\.1\\.1\\.529\\.5\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.1"
df[grep("B\\.1\\.1\\.529\\.5\\.1\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.1.x"
df[grep("B\\.1\\.1\\.529\\.5\\.1\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.1.1"
df[grep("B\\.1\\.1\\.529\\.5\\.1\\.6$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.1.6"
df[grep("B\\.1\\.1\\.529\\.5\\.1\\.23$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.1.23"
# df[grep("B\\.1\\.1\\.529\\.5\\.1\\.25$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.1.25"
# # # # df[grep("B\\.1\\.1\\.529\\.5\\.2$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.2"
# # # # df[grep("B\\.1\\.1\\.529\\.5\\.2\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.2.x"
df[grep("B\\.1\\.1\\.529\\.5\\.2\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.2.1"
df[grep("B\\.1\\.1\\.529\\.5\\.2\\.9$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.2.9"
df[grep("B\\.1\\.1\\.529\\.5\\.2\\.1\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.2.1.x (BF.x)"
# df[grep("B\\.1\\.1\\.529\\.5\\.3$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.3"
# df[grep("B\\.1\\.1\\.529\\.5\\.3\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.3.x"
# df[grep("B\\.1\\.1\\.529\\.5\\.3\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.3.1"
df[grep("B\\.1\\.1\\.529\\.5\\.3\\.1\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.3.1.x (BE.x)"
# df[grep("B\\.1\\.1\\.529\\.5\\.3\\.1\\.1\\.1\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BE.1.1.1" # There are very few of these
# df[grep("B\\.1\\.1\\.529\\.5\\.3\\.1\\.1\\.1\\.1\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BE.1.1.1.x (BQ.x)" # There are very few of these
df[grep("B\\.1\\.1\\.529\\.5\\.3\\.1\\.1\\.1\\.1\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BE.1.1.1.1 (BQ.1)"
df[grep("B\\.1\\.1\\.529\\.5\\.3\\.1\\.1\\.1\\.1\\.1\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BE.1.1.1.x (BQ.1.x)"
df[grep("B\\.1\\.1\\.529\\.5\\.3\\.1\\.1\\.1\\.1\\.1\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BE.1.1.1.1 (BQ.1.1)"
# df[grep("B\\.1\\.1\\.529\\.5\\.6$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.6"
# df[grep("B\\.1\\.1\\.529\\.5\\.6\\.",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.6.x"
df[grep("B\\.1\\.1\\.529\\.5\\.6\\.2\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.6.2.1 (BW.1)"
df[grep("B\\.1\\.1\\.529\\.5\\.6\\.2\\.1\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "BA.5.6.2.1.1 (BW.1.1)"
df[df[,"Omicron_three"]=="Omicron_Recomb","Omicron_three"] <- "Other Omicron"
df[grep("^XBB",df[,"New_Lineage"]), "Omicron_three"] <- "XBB/XBB.x"
df[grep("^XBB\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "XBB.1"
df[grep("^XBB\\.1\\.5$",df[,"New_Lineage"]), "Omicron_three"] <- "XBB.1.5"
# df[grep("^XBB$",df[,"New_Lineage"]), "Omicron_three"] <- "XBB.x" # These should not exist
# df[grep("^XBB\\.",df[,"New_Lineage"]), "Omicron_three"] <- "XBB.x"
# df[grep("^XBB\\.1$",df[,"New_Lineage"]), "Omicron_three"] <- "XBB.1"
# df[grep("^XBB\\.1\\.5$",df[,"New_Lineage"]), "Omicron_three"] <- "XBB.1.5"
df[df[,"Omicron_three"]=="BA.4.x","Omicron_three"] <- "BA.4/BA.4.x"
df[df[,"Omicron_three"]=="BA.5.x","Omicron_three"] <- "BA.5/BA.5.x"
# sort(table(df[,"Omicron_three"]))
#         Other Omicron             XBB/XBB.x     BE.1.1.1.1 (BQ.1)
#                   148                   174                   340
#     BA.5.6.2.1 (BW.1) BA.5.6.2.1.1 (BW.1.1)              BA.5.1.1
#                   426                   478                   574
#               XBB.1.5                 XBB.1             BA.5.1.23
#                   595                   607                   650
#     BA.5.3.1.x (BE.x)                BA.2.x     BA.5.2.1.x (BF.x)
#                   671                   858                   874
#   BE.1.1.1.1 (BQ.1.1)              BA.5.1.x                BA.1.x
#                   964                   995                  1079
#              BA.5.2.9                BA.2.9   BE.1.1.1.x (BQ.1.x)
#                  1133                  1237                  1350
#              BA.5.2.1                  BA.2              BA.5.1.6
#                  1412                  1424                  1492
#           BA.4/BA.4.x           BA.5/BA.5.x                  BA.1
#                  1592                  1600                  2054
#               BA.1.15                BA.5.1             BA.2.12.1
#                  2505                  2669                  3052
#                BA.1.1                Others
#                  9007                 46508


# sort(table(df[df[,"Omicron_three"]=="Other Omicron","Linaje Pangolin"])) # Get the subgroup of pango variants (change search term accordingly)
#       XAJ       XBD        XZ       XAF       XAG       XAZ       XBF       XAH
#         1         1         1         2         2         2         2         4
#       XAP       XAM B.1.1.529       XAS
#         4        12        39        78

rbd_lvl <- read.table("01_data/2023-03-21_RBD_full_lineages.tsv", header=T, quote="", sep="\t", check.names=FALSE, fill=TRUE) # This is a table of RBD mutations containing how much have accumulated per variant and gives them a level (available at https://docs.google.com/spreadsheets/d/1OTWogpyvWNTlK0ww7TlDcI4l_SkZt3Z8nYiSp2YARys/edit?pli=1#gid=0)
rbd_mut <- rbd_lvl[,"chili pepper"] # save the actual level.
temp <- names(rbd_mut) <- rbd_lvl[,"Full lineage"] # use the full lineage name for all first
temp[grep("^X",rbd_lvl[,"pango Lineage"])] <- rbd_lvl[grepl("^X",rbd_lvl[,"pango Lineage"]), "pango Lineage"] # and use the original name for recombinants only (two more are missing as they have a ">" in their names but should be avoided altogether
names(rbd_mut) <- temp # Now rename the vector (for using it as hash)
rbd_mut <- rbd_mut[names(rbd_mut)!=""] # and remove empty items
rbd_mut <- as.factor(rbd_mut) # cast as factor
# levels(rbd_mut)
# [1] ""                 "LEVEL |"          "LEVEL ||"         "LEVEL |||"
# [5] "LEVEL ||| |"      "LEVEL ||| ||"     "LEVEL ||| |||"    "LEVEL ||| ||| |"
# [9] "LEVEL ||| ||| ||"
levels(rbd_mut) <- paste0("L",0:8)
df[,"Omicron_four"] <- rbd_mut[df[,"New_Lineage"]]
df[,"Omicron_four"] <- as.character(df[,"Omicron_four"])
df[is.na(df[,"Omicron_four"]),"Omicron_four"] <- df[is.na(df[,"Omicron_four"]),"Omicron_two"]
df[df[,"Omicron_four"]=="BA.3.x","Omicron_four"] <- "Other_Omicron"
df[df[,"Omicron_four"]=="Omicron_Recomb","Omicron_four"] <- "Other_Omicron"
df[df[,"Omicron_four"]=="L0","Omicron_four"] <- "Other_Omicron"
# table(df[,"Omicron_four"])
#        BA.1.x        BA.2.x        BA.4.x        BA.5.x            L1
#          3623          4659           660          3429          2054
#            L2            L3            L4            L5            L6
#          9007          8709           572          2314          1732
#            L7 Other_Omicron        Others
#          1433          1768         46508

df[,"VOCs"] <- df[,"Delta"] # make one more column for general VOCs and VOIs (start with this one as it has the most items
df[grep("^B\\.1\\.1\\.222$",df[,"New_Lineage"]), "VOCs"] <- "222"
df[grep("^B\\.1\\.1\\.519$",df[,"New_Lineage"]), "VOCs"] <- "519"
df[grep("^B\\.1\\.1\\.7$",df[,"New_Lineage"]), "VOCs"] <- "Alpha"
df[grep("^B\\.1\\.1\\.7\\.",df[,"New_Lineage"]), "VOCs"] <- "Alpha"
df[grep("^B\\.1\\.1\\.28$",df[,"New_Lineage"]), "VOCs"] <- "Gamma"
df[grep("^B\\.1\\.1\\.28\\.",df[,"New_Lineage"]), "VOCs"] <- "Gamma"
df[df[,"Recomb"]=="Recomb","VOCs"] <- "Recomb"
df[df[,"Omicron"]=="Omicron","VOCs"] <- "Omicron"
df[grepl("^X",df[,"New_Lineage"]),"Recomb_two"] <- df[grepl("^X",df[,"New_Lineage"]),"New_Lineage"] # UPDAT 2023-03-21: We now split recombinants for a specific case
df[df[,"Recomb_two"]!="Others","VOCs"] <- df[df[,"Recomb_two"]!="Others","Recomb_two"]
# sort(table(df[,"VOCs"]),decreasing=TRUE)
#   Omicron     Delta       519    Others     Gamma     Alpha       222     XBB.1
#     38475     25706      8038      6431      2822      1853      1332       607
#   XBB.1.5        XB       XAS   XBB.1.2     XBB.2   XBB.1.4       XBB       XAM
#       595       326        78        58        43        22        15        12
#   XBB.1.6   XBB.1.9       XAH       XAP   XBB.1.1   XBB.2.2     XBB.3       XAF
#        12        11         4         4         4         4         4         2
#       XAG       XAZ       XBF       XAJ XBB.1.4.1       XBD        XZ
#         2         2         2         1         1         1         1

# sort(table(df[df[,"VOCs"]=="Others","Linaje Pangolin"])) # The following are collated into the "Others" group
#          A    A.2.5.1        A.3   B.1.1.10  B.1.1.122  B.1.1.128  B.1.1.228
#          1          1          1          1          1          1          1
#  B.1.1.231   B.1.1.33  B.1.1.334  B.1.1.367  B.1.1.403  B.1.1.416  B.1.1.434
#          1          1          1          1          1          1          1
#    B.1.1.8   B.1.1.93    B.1.147 B.1.177.53    B.1.201    B.1.206    B.1.238
#          1          1          1          1          1          1          1
#    B.1.240    B.1.258    B.1.280    B.1.333  B.1.36.31    B.1.361    B.1.409
#          1          1          1          1          1          1          1
#  B.1.438.1    B.1.503    B.1.556    B.1.567    B.1.572    B.1.576    B.1.588
#          1          1          1          1          1          1          1
#    B.1.612    B.1.623       BB.2     C.37.1   B.1.1.71    B.1.126    B.1.229
#          1          1          1          1          2          2          2
#  B.1.243.1    B.1.267    B.1.319    B.1.324    B.1.375    B.1.415    B.1.533
#          2          2          2          2          2          2          2
#    B.1.565    B.1.577    B.1.599    B.1.617    B.1.1.1  B.1.1.189  B.1.1.207
#          2          2          2          2          3          3          3
#  B.1.1.518    B.1.111    B.1.405    B.1.582       C.23    A.2.5.2       AZ.3
#          3          3          3          3          3          4          4
#  B.1.1.161  B.1.1.348    B.1.153    B.1.234    B.1.404    B.1.595    B.1.366
#          4          4          4          4          4          4          5
#  B.1.617.1        A.2  B.1.1.244  B.1.1.329  B.1.1.512  B.1.621.1  B.1.1.362
#          5          6          6          6          6          6          7
#    B.1.636        A.1        A.5  B.1.1.318    B.1.499    B.1.245    B.1.625
#          7          8          8          8          8          9          9
#    B.1.637  B.1.1.285  B.1.1.517    B.1.399    B.1.351  B.1.1.322    B.1.369
#          9         11         13         15         19         20         20
#  B.1.1.316    B.1.400    B.1.320  B.1.1.344    B.1.397    B.1.239    B.1.631
#         21         22         32         33         35         37         37
#    B.1.635    B.1.232    B.1.627    B.1.396  B.1.243.2    B.1.632    B.1.558
#         38         48         49         52         53         56         65
# Unassigned    B.1.526    B.1.241    B.1.551  B.1.415.1      A.2.5    B.1.561
#         69         70         72         73         87         92         96
#    B.1.634    B.1.189  B.1.1.432      B.1.2       C.37    B.1.427    B.1.429
#         97        111        123        215        220        257        259
#    B.1.609    B.1.621    B.1.243      B.1.1          B        B.1
#        324        437        482        544        572       1311



df[,"VOCs2"] <- df[,"VOCs"] # add a new column with the main omicron lineages specified
df[df[,"VOCs2"]=="Omicron","VOCs2"] <- df[df[,"VOCs2"]=="Omicron","Omicron_three"]
# sort(table(df[,"VOCs2"]),decreasing=TRUE)
#                 Delta                BA.1.1                   519
#                 25706                  9007                  8038
#                Others             BA.2.12.1                 Gamma
#                  6431                  3052                  2822
#                BA.5.1               BA.1.15                  BA.1
#                  2669                  2505                  2054
#                 Alpha           BA.5/BA.5.x           BA.4/BA.4.x
#                  1853                  1600                  1592
#              BA.5.1.6                  BA.2              BA.5.2.1
#                  1492                  1424                  1412
#   BE.1.1.1.x (BQ.1.x)                   222                BA.2.9
#                  1350                  1332                  1237
#              BA.5.2.9                BA.1.x              BA.5.1.x
#                  1133                  1079                   995
#   BE.1.1.1.1 (BQ.1.1)     BA.5.2.1.x (BF.x)                BA.2.x
#                   964                   874                   858
#     BA.5.3.1.x (BE.x)             BA.5.1.23                 XBB.1
#                   671                   650                   607
#               XBB.1.5              BA.5.1.1 BA.5.6.2.1.1 (BW.1.1)
#                   595                   574                   478
#     BA.5.6.2.1 (BW.1)     BE.1.1.1.1 (BQ.1)                    XB
#                   426                   340                   326
#                   XAS               XBB.1.2                 XBB.2
#                    78                    58                    43
#         Other Omicron               XBB.1.4                   XBB
#                    39                    22                    15
#                   XAM               XBB.1.6               XBB.1.9
#                    12                    12                    11
#                   XAH                   XAP               XBB.1.1
#                     4                     4                     4
#               XBB.2.2                 XBB.3                   XAF
#                     4                     4                     2
#                   XAG                   XAZ                   XBF
#                     2                     2                     2
#                   XAJ             XBB.1.4.1                   XBD
#                     1                     1                     1
#                    XZ
#                     1

df[,"Omicron_BWart"] <- df[,"Omicron"]
df[df[,"Omicron_BWart"] == "Omicron", "Omicron_BWart"] <- "Other Omicron"
df[grep("B\\.1\\.1\\.529\\.2$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BA.2"
df[grep("B\\.1\\.1\\.529\\.2\\.9$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BA.2.9"
df[grep("B\\.1\\.1\\.529\\.2\\.12\\.1$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BA.2.12.1"
df[grep("B\\.1\\.1\\.529\\.4\\.1$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BA.4.1"
df[grep("B\\.1\\.1\\.529\\.5\\.1$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BA.5.1"
df[grep("B\\.1\\.1\\.529\\.5\\.1\\.1$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BA.5.1.1"
df[grep("B\\.1\\.1\\.529\\.5\\.1\\.23$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BA.5.1.23"
df[grep("B\\.1\\.1\\.529\\.5\\.2$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BA.5.2"
df[grep("B\\.1\\.1\\.529\\.5\\.2\\.1$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BA.5.2.1"
df[grep("B\\.1\\.1\\.529\\.5\\.6$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BA.5.6"
df[grep("B\\.1\\.1\\.529\\.5\\.2\\.1\\.1",df[,"New_Lineage"]), "Omicron_BWart"] <- "BF.1"
df[grep("B\\.1\\.1\\.529\\.5\\.3\\.1\\.1\\.1\\.1\\.1$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BQ.1"
df[grep("B\\.1\\.1\\.529\\.5\\.3\\.1\\.1\\.1\\.1\\.1\\.",df[,"New_Lineage"]), "Omicron_BWart"] <- "BQ.1.x"
df[grep("B\\.1\\.1\\.529\\.5\\.3\\.1\\.1\\.1\\.1\\.1\\.1$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BQ.1.1"
df[grep("B\\.1\\.1\\.529\\.5\\.3\\.1\\.1\\.1\\.1\\.1\\.14$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BQ.1.14"
df[grep("B\\.1\\.1\\.529\\.5\\.6\\.2",df[,"New_Lineage"]), "Omicron_BWart"] <- "BA.5.6.2"
df[grep("B\\.1\\.1\\.529\\.5\\.6\\.2\\.1$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BW.1"
df[grep("B\\.1\\.1\\.529\\.5\\.6\\.2\\.1\\.1$",df[,"New_Lineage"]), "Omicron_BWart"] <- "BW.1.1"
df[grep("^XBB",df[,"New_Lineage"]), "Omicron_BWart"] <- "XBB.x"
# sort(table(df[,"Omicron_BWart"]))
#      BA.5.6.2        BA.5.6          BQ.1       BQ.1.14          BW.1
#            33           215           340           386           426
#        BA.5.2          BF.1        BW.1.1        BA.4.1      BA.5.1.1
#           455           473           478           519           574
#     BA.5.1.23        BQ.1.1        BQ.1.x        BA.2.9         XBB.x
#           650           964           964          1237          1376
#      BA.5.2.1          BA.2        BA.5.1     BA.2.12.1 Other Omicron
#          1412          1424          2669          3052         22313
#        Others
#         46508

# df <- df[df[,"Estado"]=="Yucatan",]
# df[,"Fecha de recoleccion"] <- as.Date(df[,"Fecha de recoleccion"])
# df <- df[df[,"Fecha de recoleccion"] >= "2022-05-22",]
# df <- df[df[,"Fecha de recoleccion"] < "2022-11-27",]
# sort(table(df[df[,"Omicron_BWart"]=="Other Omicron","Linaje Pangolin"])) # Get the subgroup of pango variants (change search term accordingly)


# Now, we'll be adding regions
edos <- read.table("01_data/EstadoRegion_v5.tsv",header=T, sep='\t', skip=0, comment.char='',fill=FALSE, check.names=FALSE, stringsAsFactors = FALSE, row.names=1)
df <- cbind(df, edos[df[,"Estado"],]); rownames(df) <- NULL # Append region groupings
# Fix dates
df[,"Fecha de recoleccion"] <-  as.Date(df[,"Fecha de recoleccion"],"%Y-%m-%d")
df[,"Fecha de envio"] <-  as.Date(df[,"Fecha de envio"],"%Y-%m-%d")
# Remove all those newer than the last 16 nov date
print("Last date in the whole table:")
test_date <- max(as.Date(df[,"Fecha de recoleccion"]));test_date
# [1] "2023-02-24"
# Fix gender
df[,"Genero"] <- sub("[Ff]emale","F",df[,"Genero"])
df[,"Genero"] <- sub("[Mm]ale","M",df[,"Genero"])
df[,"Genero"] <- sub("[Uu]nknown","u",df[,"Genero"])
# Fix age
df[,"Edad"] <- sub("[Uu]nknown","u",df[,"Edad"])
df[,"Edad"] <- sub("Sin informacion","999",df[,"Edad"])
df[as.numeric(df[,"Edad"]) >150,"Edad"] <- "u"
# as.character(df[,"Edad"]
# DEPRECATED: START
# # # # Fix identifiers
# # # df[,"ID Folio"] <- sub("Sin información","u",df[,"ID Folio"])
# # # df[,"SINAVE ID"] <- sub("Sin información","u",df[,"SINAVE ID"])
# # # df[,"SINOLAVE ID"] <- sub("Sin información","u",df[,"SINOLAVE ID"])
# DEPRECATED: END
# Fix type of patient
df[,"Estatus del paciente"] <- sub("[Aa]mbulatory","Amb",df[,"Estatus del paciente"])
df[,"Estatus del paciente"] <- sub("[Hh]ospitalized","Hosp",df[,"Estatus del paciente"])
df[,"Estatus del paciente"] <- sub("[Dd]eceased","Dec",df[,"Estatus del paciente"])
df[,"Estatus del paciente"] <- sub("[Uu]nknown","Unk",df[,"Estatus del paciente"])
# Append the month
df[,"Month"] <- format(as.Date(df[,"Fecha de recoleccion"]),"%Y-%m(%b)")
pdf("01_data/Monthly_genomes.pdf",width=14)
	barplot(las=2,table(df[,"Month"]), border=F, col="cornflowerblue", ylab="Total Genomes")
dev.off()
# Append the Year
df[,"Year"] <- format(as.Date(df[,"Fecha de recoleccion"]),"%Y")
# Append week
	# To do this, first create a calendar (week 1 starts on sunday of the first complete week in the year and reaches up to 52)
	week <- data.frame("week"=c(paste0("20W",sprintf('%0.2d', rep(1:52,each=7))),paste0("21W",sprintf('%0.2d', rep(1:52,each=7))),paste0("22W",sprintf('%0.2d', rep(1:52,each=7))),paste0("23W",sprintf('%0.2d', rep(1:52,each=7)))))
	# Now rename the rows to use them as a hash (dictionary)
	rownames(week) <- seq(as.Date("2020/01/05"), as.Date("2023/12/31"), by="day")[1:nrow(week)]
	# write.table(week, "week_calendar.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	# Now used the newly created dictionary to append the corresponding week
	df[,"Week"] <- week[as.character(df[,"Fecha de recoleccion"]),1]
# Next, append age categories
	# First create dictionary
	age1 <- data.frame("Age"=as.character(rep(paste0(sprintf('%0.2d',1:14),":",sapply(seq(0,130,10),function(x) paste0(x,"-",x+9))),each=10)), stringsAsFactors = FALSE)
	rownames(age1) <- 0:(nrow(age1)-1) # and use age1s as rownames
	age1[as.numeric(rownames(age1))>=100,1] <- "11:100+" # This was later added to collate all 100 an over
	df[,"Age_range_by10"] <- age1[as.character(df[,"Edad"]),1] # Now append it to the main dataframe
	# Repeat with smaller bins
	age2 <- data.frame("Age"=rep(paste0(sprintf('%0.2d',1:28),":",sapply(seq(0,135,5),function(x) paste0(x,"-",x+4))),each=5), stringsAsFactors = FALSE)
	rownames(age2) <- 0:(nrow(age2)-1) # and use age2s as rownames
	age2[as.numeric(rownames(age2))>=100,1] <- "21:100+" # This was later added to collate all 100 an over
	df[,"Age_range_by5"] <- age2[as.character(df[,"Edad"]),1] # Now append it to the main dataframe
	# Repeat with a new classification
	age3 <- data.frame("age3"=c(rep("1:0-12", 13), rep("2:13-25",13),rep("3:26-45",20), rep("4:46-60",15), rep("5:61-74",14), rep("6:75+",55)), stringsAsFactors = FALSE)
	rownames(age3) <- 0:(nrow(age3)-1) # and use age3s as rownames
	df[,"Age_range_manual"] <- age3[as.character(df[,"Edad"]),1] # Now append it to the main dataframe
	# Added a new age scheme matching vaccination
	age4 <- data.frame("age4"=c(rep("1:0-17", 18), rep("2:18-29",12),rep("3:30-39",10), rep("4:40-49",10), rep("5:50-59",10), rep("6:60+",100)), stringsAsFactors = FALSE)
	rownames(age4) <- 0:(nrow(age4)-1) # and use age3s as rownames
	df[,"Age_vac"] <- age4[as.character(df[,"Edad"]),1] #
pdf("01_data/2022-06-27_age5_genomes.pdf")
	barplot(las=2,table(df[,"Age_range_by10"]),border=F, col="coral1")
	barplot(las=2,table(df[,"Age_range_by5"]),border=F, col="coral1")
	barplot(las=2,table(df[,"Age_range_manual"]),border=F, col="coral1")
	barplot(las=2,table(df[,"Age_vac"]),border=F, col="coral1")
dev.off()
	fix_age <- function(vect){ # Ages may be missing. Fix all and convert them to character.
		vect[vect=="u"]=NA # change "u"s to NAs to prevent warnings
		vect <- sprintf('%0.3d',as.numeric(vect))
		vect[vect=="NA"]="u" # Revert NAs (they are now text) to "u"s
		return(vect)
	}
	df[,"Edad"] <- fix_age(df[,"Edad"])
	Date <- as.Date(df[,"Fecha de recoleccion"], format= "%Y-%m-%d") # Create a vector with all dates
# 	df[Date > as.Date("2021-02-18", "%Y-%m-%d"),"Last_Vaccinated"] <- "A:60+"
# 	df[Date > as.Date("2021-05-06", "%Y-%m-%d"),"Last_Vaccinated"] <- "B:50-59"
# 	df[Date > as.Date("2021-06-09", "%Y-%m-%d"),"Last_Vaccinated"] <- "C:40-49"
# 	df[Date > as.Date("2021-07-06", "%Y-%m-%d"),"Last_Vaccinated"] <- "D:30-39"
# 	df[Date > as.Date("2021-08-02", "%Y-%m-%d"),"Last_Vaccinated"] <- "E:18-29"
# 	df[Date > as.Date("2021-11-29", "%Y-%m-%d"),"Last_Vaccinated"] <- "F:15-17"
print("Resulting total items:")
dim(df)
# [1] 86468    58
write.table(df, "01_data/AllVariants.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# print("Output total items:")
# subset <- df[df["ID Folio"]!="",] # Filter only those having a folio ID
# write.table(subset, "01_data/Metadata_wID.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
