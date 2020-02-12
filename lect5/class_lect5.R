#' ---
#' title: "Class 5: Data vsiualization and graphics in R"
#' author: "Rachael.McVicar"
#' date: "January 22nd, 2020"
#' ---


plot(1:5, col="blue", typ="o")
args(plot)
?plot

# Section 2A.
# Read the data file "weight_chart.txt"

baby <- read.table("weight_chart.txt", header = TRUE)
?read.table
plot(baby, typ= "b", col="blue")
?plot
plot(baby, typ= "b", col="blue", pch=15, main="Fat babies! Woo!", xlab= "dayyzzz old", cex=1.5)



## Read a new file of mouse genome features

mouse <- read.table("feature_counts.txt", header=TRUE, sep="\t")
dotchart(mouse$Count, labels=mouse$Feature)
  ?dotchart
?barplot
barplot(mouse$Count, horiz=TRUE, col="purple",axisnames = TRUE, names.arg = NULL)
barplot(mouse$Count, horiz=TRUE, col="purple", names.arg = mouse$Feature
      )
barplot(mouse$Count, horiz=TRUE, col="purple", names.arg = mouse$Feature,las=1)      
old.par<- par()$mar

# Section 3A.
# Read the data file "male_female_counts.txt"
counts<- read.table("male_female_counts.txt", header=TRUE, sep= "\t")
?read.table
counts<- read.delim("male_female_counts.txt")
barplot(counts$Count, names.arg = counts$Sample, col=rainbow(10))

# Section 3B.
# Read and color by value "up_down_expression.txt "
genes<- read.table("up_down_expression.txt")
nrow(genes)

