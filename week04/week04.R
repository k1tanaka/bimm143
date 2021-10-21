#' ---
#' title: "Week 4: Data Visualization Lab"
#' author: "Kaito Tanaka (A15574793)"
#' date: "October 12, 2021"
#' ---

#Week 4 Data visualization Lab

#Install the package ggplot2
#install.packages("ggplot2")

#Any time I want to use this
#package I need to load it 
library(ggplot2)

View(cars)

#A quick base R plot- this is not ggplot 
plot(cars)

#Our first ggplot
#we need data + aes + geoms
ggplot(data=cars)+aes(x=speed, y=dist)+geom_point()

p<-ggplot(data=cars)+aes(x=speed, y=dist)+geom_point()

#Add a line geom with geom_line()
p+geom_line()

#Add a trend line close to the data
p+geom_smooth()

p+geom_smooth(method="lm")

#Read in our drug expression data
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)

#Q. How many genes are in this dataset
nrow(genes) 

#Q. What are the column names and how many are there
colnames(genes)
ncol(genes)

#Q. How many 'up' regulated genes
table(genes$State)

#What fraction of total genes is up-regulated
round( table(genes$State)/nrow(genes)*100,2 )

#Let's make a first plot attempt 
g<-ggplot(data=genes)+aes(x=Condition1, y=Condition2, col=State)+geom_point()

g 
#Add some color
g+scale_color_manual(values=c("blue","gray","red"))+
  labs(title="Gene expression changes",x="Control (no drug)")+
  theme_bw()

ggplot(mtcars) + aes(x=mpg, y=disp) + geom_point()


