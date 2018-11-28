# MyersMisc R package

# initialize package
install.packages("devtools")
library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)

create_package("MyersMisc")

document()

library(MyersMisc)
?MyCalPlot
?my.val.prob.ci.2
?MyValPlot
