library(areaplot)
library("plotrix")
### Institute COMPARISON ###
# ### Single variable comparisons (split by variant)
# ### DATE
# Several dates are missing and the two variants have two different histories. Thus, we should be careful with the empty categories. In this case, we'll first get the actual extant dates in a 2-way date vs lineage table
### FUNCTIONS ###
alllocations <- c("Baja California", "Baja California Sur", "Sonora", "Sinaloa", "Chihuahua", "Durango", "Coahuila", "Nuevo Leon", "Tamaulipas", "Zacatecas", "San Luis Potosi", "Aguascalientes", "Guanajuato", "Queretaro", "Hidalgo", "State of Mexico", "Mexico City", "Puebla", "Tlaxcala", "Morelos", "Nayarit", "Jalisco", "Colima", "Michoacan", "Veracruz", "Guerrero", "Oaxaca", "Tabasco", "Chiapas", "Campeche", "Yucatan", "Quintana Roo")
fix_age <- function(vect){ # Ages may be missing. Fix all and convert them to character.
	vect[vect=="u"]=NA # change "u"s to NAs to prevent warnings
	vect <- sprintf('%0.3d',as.numeric(vect))
	vect[vect=="NA"]="u" # Revert NAs (they are now text) to "u"s
	return(vect)
}
remove_1item_cols <- function(mat){ # Check if any column has <2 items and remove it if it does
	mat <- mat[,apply(mat,2,function(x) length(unique(x)))>1]
	return(mat)
}
summarize_lineage_age <- function(mat,target){ # This is a very specific case, enable it if df has this fields
	temp <- as.numeric(mat[mat$VOC_VOI==target,4])
	temp <- temp[!is.na(temp)];
	out <- c("N"=length(temp), "Mean"=round(mean(temp),2), "Median"=round(median(temp),2), "Std.Dev."=round(sd(temp),2),quantile(temp,seq(0,1,0.05)))
	return(out)
}
skewness <- function(vect) {
	avg <- mean(vect)
	a3 <- mean((vect-avg)^3)
	out <- a3/(sd(vect)^3)
	return(out)
}
kurtosis <- function(vect) {
	avg <- mean(vect)
	a4 <- mean((vect-avg)^4)
	out <- a4/(sd(vect)^4)-3
	return(out)
}
basic_statistics <- function(vect){ # This will get a test (summarized by column)
	out <- c("Mean"=round(mean(vect),2), "Median"=round(median(vect),2), "Std.Dev."=round(sd(vect),2),quantile(vect,seq(0,1,0.05)),"Skewness"=skewness(vect),"Kurtosis"=kurtosis(vect))
	return(out)
}
stats_per_column <- function(mat){ # This is a very specific case, enable it if df has this
	out <- apply(mat,2,basic_statistics)
	return(out)
}
twoway_table <- function(mat){ # Cross two variables from a 2 column matrix (var1, var2), return a table of var2 x var1 with n columns depending on total var2 items
	xtable <- xtabs(rep(1,nrow(mat))~mat[,1]+mat[,2], data=mat)
	return(xtable)
}
fill_missing_items <- function(mat,vector){ # Gets a table and a vector with all rows that should be present, outputs an expanded row collection with missing dates for complete calendar in that range. rownames should have date format as %Y-%m-%d
	xtable <- as.data.frame(matrix(0,nrow=length(vector),ncol=ncol(mat)), stringsAsFactors = FALSE) # Create empty vessel for output
	rownames(xtable) <- vector # use the vector as rownames
	colnames(xtable) <- colnames(mat) # and inherit the names
	invisible(sapply(rownames(mat),function(x) {xtable[x,] <<- mat[x,]}))	# append the original values in the corresponding places (write to higher env variable
	return(xtable)
}
rolling_N_avg <- function(mat,int){ # Input should be a table with continuous data at rows and the desired interval for the mean. If an even number is provided, the average will be placed one position to the right of as there is no single middle number
	before <- after <- trunc(int/2) # initialize with same range before and after each position
	if(int%%2==0){after <- after-1} # If even, shift the upper half of the range by 1 position (the mean will be calculated for the next value next to the middle as there is no exact number in it)
	roll <- sapply((before+1):(nrow(mat)-after),function(y) {apply(mat,2, function(x) mean(as.numeric(x[(y-before):(y+after)])))})
	colnames(roll) <- row.names(mat)[(before+1):(nrow(mat)-after)]
	roll <- t(roll)
	return(roll)
}
# ### Define functions for crossing variables
date_vs_X <- function(mat,string,int) { # Input should have at least a "Date" column with format as %Y-%m-%d, a target vector of dates that should be present (passed as a string vector) and an integer for the days in the rolling average
	dateX <- twoway_table(mat[,c("Date",string)])
	dates <- seq(as.Date(rownames(dateX)[1]),as.Date(rownames(dateX)[nrow(dateX)]),by="day") # create the complete date range
	dateX <- fill_missing_items(dateX,dates) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	dateX <- rolling_N_avg(dateX,int) # smoothen with a N-day average
	return(dateX)
}
week_vs_X <- function(mat,string) { # Input should have at least a Week column, and a target vector of all items should be present (passed as a string vector)
	weekX <- twoway_table(mat[,c("Week",string)])
	range <- sub("\\.","W",range(as.numeric(sub("W",".",mat[,"Week"])))) # get first and last items
	allweeks <- c(paste0("20W",sprintf('%0.2d', 1:52)),paste0("21W",sprintf('%0.2d',1:52)),paste0("22W",sprintf('%0.2d',1:52)),paste0("23W",sprintf('%0.2d',1:52)))
	allweeks <- allweeks[grep(range[1],allweeks):grep(range[2],allweeks)]
	weekX <- fill_missing_items(weekX,allweeks) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	return(weekX)
}
month_vs_X <- function(mat,string) { # Input should have at least a Month column, and a target vector of all items should be present (passed as a string vector)
	monthX <- twoway_table(mat[,c("Month",string)])
	range <- range(as.Date(sub("\\(.*","-01",rownames(monthX))))
	allmonths <- format(seq(range[1],range[2],by="month"),"%Y-%m(%b)")
	monthX <- fill_missing_items(monthX,allmonths) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	return(monthX)
}
age_vs_X <- function(mat,string,int) { # Input should have at least a Age column, and a target vector of all items should be present (passed as a string vector), int gets a number for a rolling avg
	ageX <- twoway_table(mat[,c("Age",string)])
	range <- range(as.numeric(rownames(ageX)[rownames(ageX)!="u"]))
	allages <- sprintf('%0.3d',seq(range[1],range[2]));if(sum(rownames(ageX)=="u")){allages[(length(allages)+1)]="u"}
	ageX <- fill_missing_items(ageX,allages) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	ageX <- rolling_N_avg(ageX,int)
	return(ageX)
}
age5_vs_X <- function(mat,string) { # Input should have at least a Age_range_by5 column, and a target vector of all items should be present (passed as a string vector)
	ageX <- twoway_table(mat[,c("Age_range_by5",string)])
	allages <- c(paste0(sprintf('%0.2d',1:20),":",sapply(seq(0,99,5),function(x) paste0(x,"-",x+4))),"21:100+","u")
	ageX <- fill_missing_items(ageX,allages) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	return(ageX)
}
age10_vs_X <- function(mat,string) { # Input should have at least a Age_range_by10 column, and a target vector of all items should be present (passed as a string vector)
	ageX <- twoway_table(mat[,c("Age_range_by10",string)])
	allages <- c(paste0(sprintf('%0.2d',1:10),":",sapply(seq(0,99,10),function(x) paste0(x,"-",x+9))),"11:100+","u")
	ageX <- fill_missing_items(ageX,allages) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	return(ageX)
}
ageM_vs_X <- function(mat,string) { # Input should have at least a Age_range_manual column, and a target vector of all items should be present (passed as a string vector)
	ageX <- twoway_table(mat[,c("Age_range_manual",string)])
	allages <- c("1:0-12", "2:13-25", "3:26-45", "4:46-60", "5:61-74", "6:75+", "u")
	ageX <- fill_missing_items(ageX,allages) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	return(ageX)
}
datediff_vs_X <- function(mat,string,int) { # Input should have at least a date column, and a target vector of all items should be present (passed as a string vector), int gets a number for a rolling avg
  daysX <- twoway_table(mat[,c("DaysPassed",string)])
  rownames(daysX) <- sprintf('%0.3d',as.numeric(rownames(daysX)))
  range <- range(as.numeric(rownames(daysX)[rownames(daysX)!="u"]))
  allages <- sprintf('%0.3d',seq(range[1],range[2]));if(sum(rownames(daysX)=="u")){allages[(length(allages)+1)]="u"}
  daysX <- fill_missing_items(daysX,allages) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
  daysX <- rolling_N_avg(daysX,int)
  #daysX <- daysX[1:99,] # subset only 0-100
  return(daysX)
}
gender_vs_X <- function(mat,string) { # Input should have at least a Gender column, and a target vector of all items should be present (passed as a string vector)
	genderX <- twoway_table(mat[,c("Gender",string)])
	allgenders <- c("F", "M", "u")
	genderX <- fill_missing_items(genderX,allgenders) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	return(genderX)
}
location_vs_X <- function(mat,string) { # Input should have at least a Gender column, and a target vector of all items should be present (passed as a string vector)
	locationX <- twoway_table(mat[,c("State",string)])
# 	alllocations <- c("Baja California", "Baja California Sur", "Sonora", "Sinaloa", "Chihuahua", "Durango", "Coahuila", "Nuevo Leon", "Tamaulipas", "Zacatecas", "San Luis Potosi", "Aguascalientes", "Guanajuato", "Queretaro", "Hidalgo", "State of Mexico", "Mexico City", "Puebla", "Tlaxcala", "Morelos", "Nayarit", "Jalisco", "Colima", "Michoacan", "Veracruz", "Guerrero", "Oaxaca", "Tabasco", "Chiapas", "Campeche", "Yucatan", "Quintana Roo") # This was made global
	locationX <- fill_missing_items(locationX,alllocations) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	return(locationX)
}
region7_vs_X <- function(mat,string) { # Input should have at least a Gender column, and a target vector of all items should be present (passed as a string vector)
	regionX <- twoway_table(mat[,c("Region_7",string)])
	allregions <- c("CN", "CS", "NE", "NW", "S", "SE", "W")
	regionX <- fill_missing_items(regionX,allregions) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	return(regionX)
}
status_vs_X <- function(mat,string) { # Input should have at least a Gender column, and a target vector of all items should be present (passed as a string vector)
	statusX <- twoway_table(mat[,c("Status",string)])
	allstatus <- c("Amb", "Dec", "Hosp", "Unk")
	statusX <- fill_missing_items(statusX,allstatus) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	return(statusX)
}
agevac_vs_X <- function(mat,string) { # Input should have at least a Gender column, and a target vector of all items should be present (passed as a string vector)
	agevacX <- twoway_table(mat[,c("Age_vac",string)])
	allagevac <- c("1:0-17", "2:18-29", "3:30-39", "4:40-49", "5:50-59", "6:60+","u")
	agevacX <- fill_missing_items(agevacX,allagevac) # use the predicted missing days to get the whole date spectrum (adding 0s when required)
	return(agevacX)
}

# ### Define functions for specific and general-purpose plotting
define_plot_scheme <- function(){
	lty <- c("solid","dashed","dotted","dotdash","longdash","twodash")
	col <- c('chartreuse3', 'cornflowerblue', 'darkgoldenrod1', 'peachpuff3','mediumorchid2', 'turquoise3', 'wheat4','slategray2',"black","coral1","aquamarine2","blue2","violetred2","palegreen3","purple3","magenta1","limegreen","darkorange2","darkgray")
	lwd <- c(1.3,2)
	pars <- expand.grid(col = col, lty = lty, lwd = lwd, stringsAsFactors = FALSE) # This will create all
	return(pars)
}
printable_dates <- function(vector){ # Gets a vector with ordered dates (fomat must have %d at the end), get the position where each month stats or at day 15 and adjust them accordingly for plotting
	dates <- unique(sort(c(grep("01$|15$",vector),length(vector),1)))
	if(length(dates)>20){dates <- unique(sort(c(grep("01$",vector),length(vector),1)))}
	if((dates[length(dates)]-dates[length(dates)-1])<10){dates <- dates[-(length(dates)-1)]} # Fix dates showing
	if((dates[2]-dates[1])<10){dates <- dates[-2]}
	return(dates)
}
plot_dates_lines <- function(mat,scheme,name,total){ # Gets an input matrix with rownames containing continuous date. This returns a plot for the current graphical device. The second parameter should be a 3-column table containing color names at col1, linetype at col2, and lwd at col3 as produced by define_plot_scheme(). A string defines the name. Total contains a number of total items for printing in the X axis
	if(length(grep("Mexico",colnames(mat))>1)){mat=t(fill_missing_items(t(mat),alllocations))} # Fix missing States if present
	oripar <- par(no.readonly=TRUE) # save default params
	cex.X=0.53
	dates <- printable_dates(rownames(mat)) # Determine where each month starts
	par(oma = c(1, 1, 1, 4)) # This is just a creative fix to plot the legend outside
	with(scheme[1:ncol(mat),],matplot(mat[,1:ncol(mat)], type = 'l', col = col, lty = lty, lwd = lwd, main = paste0(name," by date"), ylab = NA, xlab = NA, las = 2, xaxt = 'n'))
	mtext(paste("(n = ",total," total observations)"))
	axis(1, las=2, at = dates, labels = rownames(mat)[dates], cex.axis = 0.75)
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend("right", legend = colnames(mat), xpd = TRUE, horiz = FALSE, inset = c(0,0), bty = "n", col = scheme$col[1:ncol(mat)], lty = scheme$lty[1:ncol(mat)], lwd = scheme$lwd[1:ncol(mat)], cex = cex.X)
	par(oripar) # reset parameters
	return(name)
}
plot_dates_areas <- function(mat,scheme,name,total){ # Gets an input matrix with rownames containing continuous date. This returns a plot for the current graphical device. The second parameter should be a 3-column table containing color names at col1, linetype at col2, and lwd at col3 as produced by define_plot_scheme(). A string defines the name. Total contains a number of total items for printing in the X axis
	if(length(grep("Mexico",colnames(mat))>1)){mat=t(fill_missing_items(t(mat),alllocations))} # Fix missing States if present
	oripar <- par(no.readonly=TRUE) # save default params
	cex.X <- 0.53
	dates <- printable_dates(rownames(mat)) # Determine where each month starts
	par(oma = c(1, 1, 1, 4)) # This is just a creative fix to plot the legend outside
	cols <- rev(scheme[1:ncol(mat),1])
	areaplot(mat[,ncol(mat):1],col=cols,ylim=c(0,100),border=NA, las=2, main=paste0(name," by date"), ylab=NA, xaxt='n', xlab=NA)
	mtext(paste("(n = ",total," total observations)"))
	axis(1, las=2, at = dates, labels = rownames(mat)[dates], cex.axis = 0.75)
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend("right", legend = colnames(mat), xpd = TRUE, horiz = FALSE, inset = c(0,0), bty = "n", col = rev(cols), pch=15,cex = cex.X)
	par(oripar) # reset parameters
	return(name)
}
plot_lines <- function(mat,scheme,name,total){ # Gets an input matrix with rownames containing continuous date. This returns a plot for the current graphical device. The second parameter should be a 3-column table containing color names at col1, linetype at col2, and lwd at col3 as produced by define_plot_scheme(). A string defines the name. Total contains a number of total items for printing in the X axis
	if(length(grep("Mexico",colnames(mat))>1)){mat=t(fill_missing_items(t(mat),alllocations))} # Fix missing States if present
	oripar <- par(no.readonly=TRUE) # save default params
	cex.X=0.53
	par(oma = c(1, 1, 1, 4)) # This is just a creative fix to plot the legend outside
	with(scheme[1:ncol(mat),],matplot(mat[,1:ncol(mat)], type = 'l', col = col, lty = lty, lwd = lwd, main = paste0(name," by X"), ylab = NA, xlab = NA, las = 2, xaxt = 'n'))
	mtext(paste("(n = ",total," total observations)"))
	axis(1, las=2, at = 1:nrow(mat), labels = rownames(mat), cex.axis = 0.75)
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend("right", legend = colnames(mat), xpd = TRUE, horiz = FALSE, inset = c(0,0), bty = "n", col = scheme$col[1:ncol(mat)], lty = scheme$lty[1:ncol(mat)], lwd = scheme$lwd[1:ncol(mat)], cex = cex.X)
	par(oripar) # reset parameters
	return(name)
}
plot_bars <- function(mat,scheme,name,total){ # Gets an input matrix with rownames containing continuous date. This returns a plot for the current graphical device. The second parameter should be a 3-column table containing color names at col1, linetype at col2, and lwd at col3 as produced by define_plot_scheme(). A string defines the name. Total contains a number of total items for printing in the X axis
	if(length(grep("Mexico",colnames(mat))>1)){mat=t(fill_missing_items(t(mat),alllocations))} # Fix missing States if present
	oripar <- par(no.readonly=TRUE) # save default params
	cex.X <- 0.53
	par(oma = c(1, 1, 1, 4)) # This is just a creative fix to plot the legend outside
	cols <- rev(scheme[1:ncol(mat),1])
	bp <- barplot(t(mat[,ncol(mat):1]),col=cols,border=NA, las=2, main=paste0(name," vs X"), ylab=NA, xaxt='n', xlab=NA)#,ylim=c(0,100))
	mtext(paste("(n = ",total," total observations)"))
	axis(1, las=2, at = bp, labels = rownames(mat), cex.axis = 0.75)
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend("right", legend = colnames(mat), xpd = TRUE, horiz = FALSE, inset = c(0,0), bty = "n", col = rev(cols), pch=15,cex = cex.X)
	par(oripar) # reset parameters
	return(name)
}
plot_areas <- function(mat,scheme,name,total){ # Gets an input matrix with rownames containing continuous date. This returns a plot for the current graphical device. The second parameter should be a 3-column table containing color names at col1, linetype at col2, and lwd at col3 as produced by define_plot_scheme(). A string defines the name. Total contains a number of total items for printing in the X axis
	if(length(grep("Mexico",colnames(mat))>1)){mat=t(fill_missing_items(t(mat),alllocations))} # Fix missing States if present
	oripar <- par(no.readonly=TRUE) # save default params
	cex.X <- 0.53
	par(oma = c(1, 1, 1, 4)) # This is just a creative fix to plot the legend outside
	cols <- rev(scheme[1:ncol(mat),1])
	areaplot(mat[,ncol(mat):1],col=cols,ylim=c(0,100),border=NA, las=2, main=paste0(name," vs X"), ylab=NA, xaxt='n', xlab=NA)
	mtext(paste("(n = ",total," total observations)"))
	axis(1, las=2, at = 1:nrow(mat), labels = rownames(mat), cex.axis = 0.75)
	par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
	plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
	legend("right", legend = colnames(mat), xpd = TRUE, horiz = FALSE, inset = c(0,0), bty = "n", col = rev(cols), pch=15,cex = cex.X)
	par(oripar) # reset parameters
	return(name)
}
vector_to_list_by_category <- function(cat_vect, obs_vect){ # gets two vectors, one for the desired categories, one for the observed items. They should have the same number of elements.
  if(length(cat_vect)!=length(obs_vect)){return("Vectors have different total items. Aborting...")}
  # obs_vect[is.na(obs_vect)] <- "u"
  outlist <- lapply(split(obs_vect,as.factor(cat_vect)),sort)
  return(data.frame(lapply(outlist, "length<-", max(lengths(outlist)))))
}

#################################################################################################################################################################################################################################################################################################################################
### MAIN ###
# ### Whole table ###
df <- read.table("01_data/AllVariants.tsv",header=T, sep='\t', skip=0, comment.char="",fill=T, check.names=FALSE, stringsAsFactors = FALSE) # columns regarding mutations were removed as some EOF error arose
# names(df)[which(names(df)=="Fecha de recoleccion")] <- "Date" # In case only some labels need updating
dim(df)
# [1] 84750    54
df <- df[,c("Accession ID", "Nombre", "Fecha de recoleccion", "Fecha de envio", "Estado", "Municipio", "Genero", "Edad","E","RdRP","Estatus del paciente","Clado Nexstrain", "New_Lineage", "Delta", "Delta_sub", "Omicron","Omicron_sub", "Omicron_two", "Omicron_three", "Omicron_four", "Omicron_BWart", "VOCs", "VOCs2", "Region_7", "Region_5", "Region_4", "Region_3", "1999_Inegi", "CEEY", "Regiones_operativas", "Om2_A", "Om2_B", "Om2_C", "Month", "Year", "Week", "Age_range_by10", "Age_range_by5","Age_range_manual", "Age_vac")]
dim(df)
# [1] 76689    43
# write.table(cbind("Old_name"=colnames(df),"New_name"=colnames(df)), "01_data/dummyCols_forAnalysis.tsv", sep="\t", col.names=TRUE, quote = FALSE, row.names=FALSE)
rename <- read.table("01_data/rename_fields_v3.tsv",header=T, sep='\t', skip=0, comment.char='',fill=F, check.names=FALSE, stringsAsFactors = FALSE, row.names=1)
names(df) <- rename[names(df),]
df[is.na(df)]="u" # Fill missing items
# Remove known bad registries
# df <- df[!(df[,"Delta_or_not"]=="Delta" & (as.Date(df[,"Date"]) < as.Date("2021-04-21", "%Y-%m-%d"))),] # Delta samples before May 2020
df <- remove_1item_cols(df) # if any column doesn't have at least 2 items, remove it (first "official" date was Apr 21th, 2021)
df[,grep("Age$", colnames(df))] <- fix_age(df[,grep("Age$", colnames(df))]) # Transform u in ages to NAs (by coercion while casting to numeric)
df[df=="u"]=NA # change "u"s to NAs to prevent warnings # UNCOMMENT if required
df[df=="Unk"]=NA
df[grep("IBT",df[,"Genome"]),"CoViGen"] <- "CoViGen"
df[grep("LANGEBIO",df[,"Genome"]),"CoViGen"] <- "CoViGen"
df[grep("CIAD",df[,"Genome"]),"CoViGen"] <- "CoViGen"
# df[grep("INER",df[,"Genome"]),"CoViGen"] <- "CoViGen"
df[grep("INER.IMSS",df[,"Genome"]),"CoViGen"] <- "CoViGen"
# df[grep("InDRE",df[,"Genome"]),"CoViGen"] <- "Others"
# df[grep("INER.*IBT",df[,"Genome"]),"CoViGen"] <- "Others"
# df[grep("INER.*INC",df[,"Genome"]),"CoViGen"] <- "Others"
# df[grep("INMEGEN",df[,"Genome"]),"CoViGen"] <- "Others"
df[is.na(df["CoViGen"]),"CoViGen"] <- "Others"

df[grep("InDRE",df[,"Genome"]),"Instituto"] <- "InDRE"
df[grep("IBT",df[,"Genome"]),"Instituto"] <- "zIBT"
df[grep("INER.IMSS",df[,"Genome"]),"Instituto"] <- "zINER.IMSS"
df[grep("LANGEBIO",df[,"Genome"]),"Instituto"] <- "zLANGEBIO"
df[grep("CIAD",df[,"Genome"]),"Instituto"] <- "zCIAD"
df[grep("ALSR",df[,"Genome"]),"Instituto"] <- "ALSR"
#df[grep("ISSEA",df[,"Genome"]),"Instituto"] <- "ISSEA" #<100 when filtered by those <100 days passed for gisaid upload
df[grep("MEXBV",df[,"Genome"]),"Instituto"] <- "MEXBV"
#df[grep("NHRC",df[,"Genome"]),"Instituto"] <- "NHRC"
df[grep("SEARCH",df[,"Genome"]),"Instituto"] <- "ALSR"
df[grep("CABANA",df[,"Genome"]),"Instituto"] <- "CABANA"
df[grep("CEMANAV",df[,"Genome"]),"Instituto"] <- "CEMANAV"
#df[grep("HIMFG",df[,"Genome"]),"Instituto"] <- "HIMFG"
df[grep("INMEGEN",df[,"Genome"]),"Instituto"] <- "INMEGEN"
df[grep("LESPNL",df[,"Genome"]),"Instituto"] <- "LESPNL"
#df[grep("LABOPAT",df[,"Genome"]),"Instituto"] <- "LABOPAT"
#df[grep("CICESE",df[,"Genome"]),"Instituto"] <- "CICESE"
#df[grep("UASLP",df[,"Genome"]),"Instituto"] <- "UASLP"
df[grep("UANL",df[,"Genome"]),"Instituto"] <- "UANL"
#df[grep("LANI",df[,"Genome"]),"Instituto"] <- "LANIIA"
df[is.na(df["Instituto"]),"Instituto"] <- "zzOthers"
df[,"DaysPassed"] <- as.Date(df[,"Date2"])-as.Date(df[,"Date"])

#df <- df[df[,"CoViGen"]=="CoViGen",]
df <- remove_1item_cols(df) # if any column doesn't have at least 2 items, remove it

# table(df["Instituto"])

dfresp <- df
dfresp[,"DaysPassed"] <- as.numeric(dfresp[,"DaysPassed"])
# dfresp[dfresp[,"DaysPassed"]>500,]
Ins <- split(dfresp,dfresp[,"Instituto"])
Ins <- lapply(Ins, function(x){vector_to_list_by_category(x[,"Week"], x[,"DaysPassed"])})
Ins <- lapply(Ins,colMeans, na.rm=T)
allnames <- unique(sort(unlist(lapply(Ins,names))))
Institutos <- names(Ins)
outmat <- matrix(0,nrow = length(Institutos), ncol=length(allnames));colnames(outmat) <- allnames; rownames(outmat) <- Institutos
for(inst in names(Ins)){
  print(inst)
  mat <- Ins[[inst]]
  for(time in names(mat)){
    #print(time)
    outmat[inst,time] <- mat[time]
  }
}
colnames(outmat) <- sub("^X","",colnames(outmat))
rownames(outmat) <- gsub("z","",rownames(outmat))
write.table(t(outmat), "days_institute_per_week_all.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

dfresp <- df
dfresp[,"DaysPassed"] <- as.numeric(dfresp[,"DaysPassed"])
Ins <- split(dfresp,dfresp[,"Instituto"])
Ins <- lapply(Ins, function(x){vector_to_list_by_category(x[,"Month"], x[,"DaysPassed"])})
Ins <- lapply(Ins,colMeans, na.rm=T)
allnames <- unique(sort(unlist(lapply(Ins,names))))
Institutos <- names(Ins)
outmat <- matrix(0,nrow = length(Institutos), ncol=length(allnames));colnames(outmat) <- allnames; rownames(outmat) <- Institutos
for(inst in names(Ins)){
  print(inst)
  mat <- Ins[[inst]]
  for(time in names(mat)){
    #print(time)
    outmat[inst,time] <- mat[time]
  }
}
colnames(outmat) <- sub("^X","",colnames(outmat))
rownames(outmat) <- gsub("z","",rownames(outmat))
write.table(t(outmat), "days_institute_per_month_all.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

# df <- df[df[,"DaysPassed"]<=100,]
dim(df)
# [1] 84750    41
# UPDATE 2022-09-23: Round Cts into integers to reduce categories
# df["Ct_E"] <- round(as.numeric(unlist(df["Ct_E"])))
# df["Ct_RdRp"] <- round(as.numeric(unlist(df["Ct_RdRp"])))

ignored <- c("Acc", "Genome", "ID_folio", "SINAVE", "SINOLAVE")
df <- df[,-which(names(df) %in% ignored)] # Drop items that should be ignored
ignored_for_date <- c("Date", "Date2", "Date_Batch", "Week", "Month", "Year", "Lineage","Municipality", "Last_Vaccinated","DaysPassed","Omicron_sub")
ignored_for_ages <- c("Age", "Age_range_by5", "Age_range_by10", "Age_range_manual", "Age_vac", "Date", "Date2", "Lineage","Municipality","DaysPassed","Omicron_sub")
ignored_for_days <- c("Age","Date", "Date2", "Lineage","Municipality","DaysPassed","Omicron_sub")

dir.create("Raw_tables",showWarnings=FALSE)
days <- 7
dateXall <- sapply(setdiff(colnames(df),ignored_for_date), function(x) {date_vs_X(df,x,days)}) # Exclude named columns. Create a set of tables for the whole table smoothen by a 7-day rolling average
# sapply(colnames(df[,-1]), function(x) {out <- date_vs_X(df,x,days);write.table(out, paste0("date_",x,".tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)}) # Output the tables
dateXall_rel <- lapply(dateXall,function(x) prop.table(x,2)) # calculate relative numbers by column (dynamic frequency per date)
dateXall_stacked <- lapply(dateXall,function(x) prop.table(x,1)*100) # Get % per item
dateXall_sums <- as.data.frame(rowSums(dateXall$State));colnames(dateXall_sums)="All_genomes" # add a sum for a general table
params <- define_plot_scheme()
CompDateXtotals <- dateXall_sums; # save a copy for comparison with national data
names(CompDateXtotals) <- "GS" # update the name to identify surveillance data
pdf("00_date_national_dates_sums_7daywindow.pdf",width=14)
	plot_dates_lines(dateXall_sums,params,"State",nrow(df))
dev.off()
pdf("01_date_continuous_dates_rawtab_7daywindow.pdf",width=14)
	for(i in names(dateXall)){print(i);plot_dates_lines(dateXall[[i]],params,i,sum(!is.na(df[,i])))}
dev.off()
pdf("02_date_freq_continuous_dates_rawtab_7daywindow.pdf",width=14)
	for(i in names(dateXall_rel)){print(i);plot_dates_lines(dateXall_rel[[i]],params,i,sum(!is.na(df[,i])))}
dev.off()
pdf("03_date_perc_continuous_dates_rawtab_7daywindow.pdf",width=14)
	for(i in names(dateXall_stacked)){print(i);plot_dates_areas(dateXall_stacked[[i]],params,i,sum(!is.na(df[,i])))}
dev.off()
# Output tables
for(i in names(dateXall)){write.table(dateXall[[i]], paste0("Raw_tables/01_date_all_lineages_and_dates_7dAvg_",i,".tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)}


# Pie charts
df[df[,"Instituto"]=="UANL","Instituto"] <- "zzOthers"
df[df[,"Instituto"]=="CABANA","Instituto"] <- "zzOthers"
df[df[,"Instituto"]=="LESPNL","Instituto"] <- "zzOthers"
inst <- df[,"Instituto"]
table(inst)
# ALSR      InDRE    INMEGEN      zCIAD       zIBT zINER.IMSS  zLANGEBIO   zzOthers 
# 3558      14883      25798       2011      12779       3487       9739       4434 
inst_covigen <- inst
inst_covigen[inst_covigen=="zIBT"] <- "CoViGen"
inst_covigen[inst_covigen=="zLANGEBIO"] <- "CoViGen"
inst_covigen[inst_covigen=="zINER.IMSS"] <- "CoViGen"
inst_covigen[inst_covigen=="zCIAD"] <- "CoViGen"
table(inst_covigen)
# ALSR  CoViGen    InDRE  INMEGEN zzOthers 
# 3558    28016    14883    25798     4434 
round(table(inst_covigen)/sum(table(inst_covigen))*100,2)
# ALSR  CoViGen    InDRE  INMEGEN zzOthers 
# 4.64    36.53    19.41    33.64     5.78 
temp <- sort(table(inst_covigen), decreasing=TRUE)
cols <- c("purple","firebrick","limegreen","royalblue","deeppink","coral1","aquamarine", "gold2","cornflowerblue","purple","wheat", "cyan4", "plum1", "chartreuse3", "azure4", "deeppink", "nabyblue","tan2", "palegreen", "brown2")
pdf("Pie_CoVigen_vs_others.pdf")
	pie3D(temp, col=cols, theta=1,shade=1, mar=c(4,4,4,4), explode=0.1, labels=names(temp))
dev.off()
inst_covigen[grep("^z", inst)] <- inst[grep("^z", inst)]
temp <- table(inst_covigen)
temp <- temp[c(5,7,6,4,3,2,8,1)]
# cols <- c("aquamarine","coral1","gold2","tan2","firebrick","limegreen","royalblue","deeppink","coral1","aquamarine", "gold2","cornflowerblue","purple","wheat", "cyan4", "plum1", "chartreuse3", "azure4", "deeppink", "nabyblue","tan2", "palegreen", "brown2")
cols <- c("darkorchid1","darkorchid2","darkorchid3","darkorchid4","firebrick","limegreen","royalblue","deeppink","coral1","aquamarine", "gold2","cornflowerblue","purple","wheat", "cyan4", "plum1", "chartreuse3", "azure4", "deeppink", "nabyblue","tan2", "palegreen", "brown2")
pdf("Pie_CoVigen_vs_others_all.pdf")
	pie3D(temp, col=cols, theta=1,shade=1, mar=c(4,4,4,4), explode=0.1, labels=names(temp))
dev.off()
temp
# zIBT  zLANGEBIO zINER.IMSS      zCIAD    INMEGEN      InDRE   zzOthers       ALSR 
# 14673      11184       3487       2498      27694      15918       5677       3619
round(temp/sum(temp)*100,2)
#  zIBT  zLANGEBIO zINER.IMSS      zCIAD    INMEGEN      InDRE   zzOthers       ALSR 
# 17.31      13.20       4.11       2.95      32.68      18.78       6.70       4.27


# ### STATES only
states <- dateXall$State
params <- define_plot_scheme()
printDat <- printable_dates(rownames(states))
dir.create("Estados", showWarnings=FALSE)
for(i in 1:32){
	print(i)
	nam <- colnames(states)[i]
	pdf(paste0("Estados/",nam,".pdf"),width=20)
	plot(las=2, states[,i], col=params[i,1], type='l', lwd=2, main=paste0("State: ", nam, " Total N: "), ylab="Average Daily Genomes", xaxt='n', xlab='')
	axis(1, las=2, at=printDat, labels=rownames(states)[printDat], cex.axis=0.8)
	dev.off()
}
write.table(states, "States_cases.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

states <- week_vs_X(df, "State")
states <- states[1:(nrow(states)-1),]
params <- define_plot_scheme()
dir.create("Estados_week", showWarnings=FALSE)
for(i in 1:32){
	print(i)
	nam <- colnames(states)[i]
	pdf(paste0("Estados_week/",nam,".pdf"),width=20)
	plot(las=2, states[,i], col=params[i,1], type='l', lwd=2, main=paste0("State: ", nam, " Total N: "), ylab="Weekly Genomes", xlab='', xaxt='n')
	axis(1, las=2, at=1:nrow(states), labels=rownames(states), cex.axis=0.8)
	dev.off()
}
write.table(states, "States_cases-week.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

# UPDATE 2023-02-11: Now with the state and the Omicron three set
df[,"State_Om3"] <- paste(df[,"State"],df[,"Omicron_three"],sep="|")
dateVsOm2COm3 <- date_vs_X(df,"State_Om3",7)
write.table(dateVsOm2COm3, "State_Om3_cases-day.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
weekVsOm2COm3 <- week_vs_X(df, "State_Om3")
write.table(weekVsOm2COm3, "State_Om3_cases-week.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
df[,"State_Om4"] <- paste(df[,"State"],df[,"Omicron_four"],sep="|")
dateVsOm2COm4 <- date_vs_X(df,"State_Om4",7)
write.table(dateVsOm2COm4, "State_Om4_cases-day.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
weekVsOm2COm4 <- week_vs_X(df, "State_Om4")
write.table(weekVsOm2COm4, "State_Om4_cases-week.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
df[,"State_OmABWart"] <- paste(df[,"State"],df[,"Omicron_BWart"],sep="|")
dateVsStOmBWart <- date_vs_X(df,"State_OmABWart",7)
write.table(dateVsStOmBWart, "State_OmBWart_cases-day.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
weekVsStOmBWart <- week_vs_X(df, "State_OmABWart")
write.table(weekVsStOmBWart, "State_OmBWart_cases-week.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


# days <- 30
# dateXall <- sapply(setdiff(colnames(df),ignored_for_date), function(x) {date_vs_X(df,x,days)}) # Exclude named columns. Create a set of tables for the whole table smoothen by a multi-day rolling average
# dateXall_sums <- as.data.frame(rowSums(dateXall$State));colnames(dateXall_sums)="GS" # add a sum for a general table
# CompDateXtotals <- dateXall_sums; # save a copy for comparison with national data
# CompDateXRegions <- dateXall[["Region_7"]]; colnames(CompDateXRegions) <- paste0("GS_",colnames(CompDateXRegions)) # extract some tables for comparison with national data
# CompDateXStatus <- dateXall[["Status"]]; colnames(CompDateXStatus) <- paste0("GS_",colnames(CompDateXStatus))
# CompDateXAge <- dateXall[["Age_vac"]]; colnames(CompDateXAge) <- paste0("GS_",colnames(CompDateXAge))

# Now in ages
df <- df[,1:38] # Remove extra columns added for the last part
ages <- 5
ignored_for_ages
ageXall <- sapply(setdiff(colnames(df),ignored_for_ages), function(x) {age_vs_X(df,x,ages)})
ageXall_rel <- lapply(ageXall,function(x) prop.table(x,2)) # calculate relative numbers by column (dynamic frequency per date)
ageXall_stacked <- lapply(ageXall,function(x) prop.table(as.matrix(x),1)*100) # Get % per item
pdf("01_age_continuous_dates_rawtab_5ywindow.pdf",width=14)
	for(i in names(ageXall)){print(i);plot_lines(ageXall[[i]],params,i,sum(!is.na(df[,i])))}
dev.off()
pdf("02_age_freq_continuous_dates_rawtab_5daywindow.pdf",width=14)
	for(i in names(ageXall_rel)){print(i);plot_lines(ageXall_rel[[i]],params,i,sum(!is.na(df[,i])))}
dev.off()
ageXall_sums <- as.data.frame(rowSums(ageXall$State));colnames(ageXall_sums)="All_genomes" # add a sum for a general table

pdf("00_age_national_ages_sums_5ywindow.pdf",width=14)
	plot_lines(ageXall_sums,params,"State",nrow(df))
dev.off()
# Output tables
for(i in names(ageXall)){write.table(ageXall[[i]], paste0("Raw_tables/01_raw_all_lineages_and_dates_30dAvg_",i,".tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)}

days <- 1
dateXall <- sapply(setdiff(colnames(df),ignored_for_date), function(x) {date_vs_X(df,x,days)})
for(i in names(dateXall)){write.table(dateXall[[i]], paste0("Raw_tables/00_date_all_raw_",i,".tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)}
ages <- 1
ageXall <- sapply(setdiff(colnames(df),ignored_for_ages), function(x) {age_vs_X(df,x,ages)})
for(i in names(ageXall)){write.table(ageXall[[i]], paste0("Raw_tables/00_age_all_raw_",i,".tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)}
ageXall_sums <- as.data.frame(rowSums(ageXall$State));colnames(ageXall_sums)="All_genomes" # add a sum for a general table
pdf("00_age_national_ages_sums_1ywindow.pdf",width=14)
  plot_bars(ageXall_sums,params,"State",nrow(df))
dev.off()
# CompAgeXtotals <- ageXall_sums # save a copy for comparison with national data
# names(CompAgeXtotals) <- "GS" # update the name to identify surveillance data
# CompAgeXRegions <- ageXall[["Region_7"]]; colnames(CompAgeXRegions) <- paste0("GS_",colnames(CompAgeXRegions)) # extract some tables for comparison with national data
# CompAgeXStatus <- ageXall[["Status"]]; colnames(CompAgeXStatus) <- paste0("GS_",colnames(CompAgeXStatus))
# CompAgeXWeeks <- ageXall[["Week"]]; colnames(CompAgeXWeeks) <- paste0("GS_",colnames(CompAgeXWeeks))
# CompAgeXMonths <- ageXall[["Month"]]; colnames(CompAgeXMonths) <- paste0("GS_",colnames(CompAgeXMonths))
# original_date_range <- range(as.Date(df$Date))
# # And add one more for the months
# MonthXall <- lapply(sapply(setdiff(colnames(df),ignored_for_date), function(x) {month_vs_X(df,x)}),function(x) x)
# # monthXall <- lapply(sapply(setdiff(colnames(df),ignored_for_date), function(x) {month_vs_X(df,x)}),function(x) prop.table(as.matrix(x),1)*100)
# MonthXall_sums <- as.data.frame(rowSums(MonthXall$State));colnames(MonthXall_sums)="GS" # add a sum for a general table
# CompMonthXtotals <- MonthXall_sums
# CompMonthXRegions <- MonthXall[["Region_7"]]; colnames(CompMonthXRegions) <- paste0("GS_",colnames(CompMonthXRegions))
# CompMonthXStatus <- MonthXall[["Status"]]; colnames(CompMonthXStatus) <- paste0("GS_",colnames(CompMonthXStatus))
# CompMonthXAge <- MonthXall[["Age_vac"]]; colnames(CompMonthXAge) <- paste0("GS_",colnames(CompMonthXAge))

# # Now in days passed
# days <- 5
# daysXall <- sapply(setdiff(colnames(df),ignored_for_days), function(x) {datediff_vs_X(df,x,days)})
# daysXall_stacked <- lapply(daysXall,function(x) prop.table(as.matrix(x),1)*100) # Get % per item
# pdf("01_days_continuous_days_rawtab_5ywindow.pdf",width=14)
# for(i in names(daysXall)){print(i);plot_lines(daysXall[[i]],params,i,sum(!is.na(df[,i])))}
# dev.off()
# 
# 
# dfresp <- df
# dfresp[,"DaysPassed"] <- as.numeric(dfresp[,"DaysPassed"])
# Ins <- split(dfresp,dfresp[,"Instituto"])
# Ins <- lapply(Ins, function(x){vector_to_list_by_category(x[,"Week"], x[,"DaysPassed"])})
# Ins <- lapply(Ins,colMeans, na.rm=T)
# allnames <- unique(sort(unlist(lapply(Ins,names))))
# Institutos <- names(Ins)
# outmat <- matrix(0,nrow = length(Institutos), ncol=length(allnames));colnames(outmat) <- allnames; rownames(outmat) <- Institutos
# for(inst in names(Ins)){
#   print(inst)
#   mat <- Ins[[inst]]
#   for(time in names(mat)){
#     #print(time)
#     outmat[inst,time] <- mat[time]
#   }
# }
# colnames(outmat) <- sub("^X","",colnames(outmat))
# rownames(outmat) <- gsub("z","",rownames(outmat))
# write.table(t(outmat), "days_institute_per_week.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
# 
# dfresp <- df
# dfresp[,"DaysPassed"] <- as.numeric(dfresp[,"DaysPassed"])
# Ins <- split(dfresp,dfresp[,"Instituto"])
# Ins <- lapply(Ins, function(x){vector_to_list_by_category(x[,"Month"], x[,"DaysPassed"])})
# Ins <- lapply(Ins,colMeans, na.rm=T)
# allnames <- unique(sort(unlist(lapply(Ins,names))))
# Institutos <- names(Ins)
# outmat <- matrix(0,nrow = length(Institutos), ncol=length(allnames));colnames(outmat) <- allnames; rownames(outmat) <- Institutos
# for(inst in names(Ins)){
#   print(inst)
#   mat <- Ins[[inst]]
#   for(time in names(mat)){
#     #print(time)
#     outmat[inst,time] <- mat[time]
#   }
# }
# colnames(outmat) <- sub("^X","",colnames(outmat))
# rownames(outmat) <- gsub("z","",rownames(outmat))
# write.table(t(outmat), "days_institute_per_month.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

# #again, without City and State of Mexico
# 
# dfresp <- df
# df <- df[-grep("Mexico",df[,"State"]),]
# df <- remove_1item_cols(df) # if any column doesn't have at least 2 items, remove it (first "official" date was Apr 21th, 2021)
# dfresp[,"DaysPassed"] <- as.numeric(dfresp[,"DaysPassed"])
# Ins <- split(dfresp,dfresp[,"Instituto"])
# Ins <- lapply(Ins, function(x){vector_to_list_by_category(x[,"Week"], x[,"DaysPassed"])})
# Ins <- lapply(Ins,colMeans, na.rm=T)
# allnames <- unique(sort(unlist(lapply(Ins,names))))
# Institutos <- names(Ins)
# outmat <- matrix(0,nrow = length(Institutos), ncol=length(allnames));colnames(outmat) <- allnames; rownames(outmat) <- Institutos
# for(inst in names(Ins)){
#   print(inst)
#   mat <- Ins[[inst]]
#   for(time in names(mat)){
#     #print(time)
#     outmat[inst,time] <- mat[time]
#   }
# }
# colnames(outmat) <- sub("^X","",colnames(outmat))
# rownames(outmat) <- gsub("z","",rownames(outmat))
# write.table(t(outmat), "days_institute_per_week_except_CDMX_and_EdoMex.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
# 
# dfresp <- df
# dfresp[,"DaysPassed"] <- as.numeric(dfresp[,"DaysPassed"])
# Ins <- split(dfresp,dfresp[,"Instituto"])
# Ins <- lapply(Ins, function(x){vector_to_list_by_category(x[,"Month"], x[,"DaysPassed"])})
# Ins <- lapply(Ins,colMeans, na.rm=T)
# allnames <- unique(sort(unlist(lapply(Ins,names))))
# Institutos <- names(Ins)
# outmat <- matrix(0,nrow = length(Institutos), ncol=length(allnames));colnames(outmat) <- allnames; rownames(outmat) <- Institutos
# for(inst in names(Ins)){
#   print(inst)
#   mat <- Ins[[inst]]
#   for(time in names(mat)){
#     #print(time)
#     outmat[inst,time] <- mat[time]
#   }
# }
# colnames(outmat) <- sub("^X","",colnames(outmat))
# rownames(outmat) <- gsub("z","",rownames(outmat))
# write.table(t(outmat), "days_institute_per_month_except_CDMX_and_EdoMex.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
# 
# 
# 
# ###########################################################################################################
# # ### CoViGen ONLY
# subtime <- df # Start with a copy (if we ever need to get back to the original table)
# subtime[subtime=="u"]=NA # change "u"s to NAs to prevent warnings # UNCOMMENT if required
# # subtime[subtime=="Unk"]=NA
# subtime$Date <- as.Date(subtime$Date, format= "%Y-%m-%d") # Make sure casting to Date type was done
# # subtime <- subset(subtime, Date >= "2021-12-01" & Date <= "2022-04-15")
# # subtime <- subtime[subtime["Omicron_or_not"]=="Omicron",]
# subtime <- subtime[subtime[,"CoViGen"]=="CoViGen",]
# subtime <- remove_1item_cols(subtime) # if any column doesn't have at least 2 items, remove it
# dim(subtime)
# # [1] 21947    39
# days <- 7
# dateXsubtime <- sapply(setdiff(colnames(subtime),ignored_for_date), function(x) {date_vs_X(subtime,x,days)})
# dateXsubtime_rel <- lapply(dateXsubtime,function(x) prop.table(x,2)) # calculate relative numbers by column (dinamic frequency per date)
# dateXsubtime_stacked <- lapply(dateXsubtime,function(x) prop.table(x,1)*100) # Get % per item
# dateXsubtime_sums <- as.data.frame(rowSums(dateXsubtime$State));colnames(dateXsubtime_sums)="All_genomes" # add a sum for a general table
# params <- define_plot_scheme()
# pdf("00_date_national_dates_sums_dec-apr-7daywindow.pdf",width=14)
# plot_dates_lines(dateXsubtime_sums,params,"State",nrow(subtime))
# dev.off()
# params <- define_plot_scheme()
# pdf("04_date_continuous_dates_dec-apr_7daywindow.pdf",width=14)
# for(i in names(dateXsubtime)){print(i);plot_dates_lines(dateXsubtime[[i]],params,i,sum((!is.na(subtime[,i]))*(!is.na(subtime[,"Date"]))))}
# dev.off()
# pdf("05_date_freq_continuous_dates_dec-apr_7daywindow.pdf",width=14)
# for(i in names(dateXsubtime_rel)){print(i);plot_dates_lines(dateXsubtime_rel[[i]],params,i,sum((!is.na(subtime[,i]))*(!is.na(subtime[,"Date"]))))}
# dev.off()
# pdf("06_date_perc_continuous_dates_dec-apr_7daywindow.pdf",width=14)
# for(i in names(dateXsubtime_stacked)){print(i);plot_dates_areas(dateXsubtime_stacked[[i]],params,i,sum((!is.na(subtime[,i]))*(!is.na(subtime[,"Date"]))))}
# dev.off()
# # Output tables
# for(i in names(dateXsubtime)){write.table(dateXsubtime[[i]], paste0("Raw_tables/04_date_all_lineages-dec-apr_7dAvg_",i,".tsv"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)}
# # # Again, for comparison with the complete surveillance table
# # days <- 30
# # dateXsubtime <- sapply(setdiff(colnames(subtime),ignored_for_date), function(x) {date_vs_X(subtime,x,days)})
# # dateXsubtime_sums <- as.data.frame(rowSums(dateXsubtime$State));colnames(dateXsubtime_sums)="GS_Delta" # add a sum for a general table
# # # append to the comparison tables
# # CompDateXtotals <- merge(CompDateXtotals,dateXsubtime_sums, by="row.names", all=TRUE);CompDateXtotals[is.na(CompDateXtotals)] <- 0; rownames(CompDateXtotals) <- CompDateXtotals[,1];CompDateXtotals <- CompDateXtotals[,-1]
# # temp <- dateXsubtime[["Region_7"]]; colnames(temp) <- paste0("GS_Delta_",colnames(temp)); CompDateXRegions <- merge(CompDateXRegions,temp, by="row.names", all=TRUE);CompDateXRegions[is.na(CompDateXRegions)] <- 0; rownames(CompDateXRegions) <- CompDateXRegions[,1];CompDateXRegions <- CompDateXRegions[,-1]
# # temp <- dateXsubtime[["Status"]]; colnames(temp) <- paste0("GS_Delta_",colnames(temp)); CompDateXStatus <- merge(CompDateXStatus,temp, by="row.names", all=TRUE);CompDateXStatus[is.na(CompDateXStatus)] <- 0; rownames(CompDateXStatus) <- CompDateXStatus[,1];CompDateXStatus <- CompDateXStatus[,-1]
# # temp <- dateXsubtime[["Age_vac"]]; colnames(temp) <- paste0("GS_Delta_",colnames(temp)); CompDateXAge <- merge(CompDateXAge,temp, by="row.names", all=TRUE);CompDateXAge[is.na(CompDateXAge)] <- 0; rownames(CompDateXAge) <- CompDateXAge[,1];CompDateXAge <- CompDateXAge[,-1]
# 
# # Now in ages
# ages <- 5
# ageXsubtime <- sapply(setdiff(colnames(subtime),ignored_for_ages), function(x) {age_vs_X(subtime,x,ages)})
# ageXsubtime_stacked <- lapply(ageXsubtime,function(x) prop.table(as.matrix(x),1)*100) # Get % per item
# ageXsubtime_sums <- as.data.frame(rowSums(ageXsubtime$State));colnames(ageXsubtime_sums)="All_genomes" # add a sum for a general table
# pdf("00_age_national_ages_sums_dec-apr-5ywindow.pdf",width=14)
# plot_lines(ageXsubtime_sums,params,"State",sum(!is.na(subtime[,"Age"])))
# dev.off()
# pdf("04_age_continuous_ages_dec-apr_5ywindow.pdf",width=14)
# for(i in names(ageXsubtime)){print(i);plot_lines(ageXsubtime[[i]],params,i,sum((!is.na(subtime[,i]))*(!is.na(subtime[,"Age"]))))}
# dev.off()
# ages <- 5
# temp <- as.numeric(subtime[,"Age"])
# temp[is.na(temp)]=10000
# age_lt90 <- subtime[temp<90,]
# ageXsubtime_lt90 <- sapply(setdiff(colnames(age_lt90),ignored_for_ages), function(x) {age_vs_X(age_lt90,x,ages)})
# ageXsubtime_stacked_lt90 <- lapply(ageXsubtime_lt90,function(x) prop.table(as.matrix(x),1)*100) # Get % per item
# pdf("06_age_continuous_ages_dec-apr_5ywindow_perc_lt90.pdf",width=14)
# for(i in names(ageXsubtime_stacked_lt90)){print(i);plot_lines(ageXsubtime_stacked_lt90[[i]],params,i,sum((!is.na(age_lt90[,i]))*(!is.na(age_lt90[,"Age"]))))}
# dev.off()
# save.image("checkpoint1.Rdata")
