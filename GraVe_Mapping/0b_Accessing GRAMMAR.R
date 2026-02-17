# Way 1 - doesn't work - not on CRAN
#install.packages("GenABEL",dependencies=c("Depends","Suggests"))
#av <- available.packages(filters=list())
#av[av[, "Package"] == 'GenABEL', ]

# Way 2 - manual installation - doesn't work
#devtools::install("C:\\Users\\jtanshengyi\\Downloads\\GenABEL_1.8-0\\GenABEL")

# Way 3 - combining Rtools and devtools - seems to work!
# devtools
#install.packages('devtools')
library(devtools)

# Rtools
#install.packages('installr')
library(installr)
#install.Rtools()

# Make the files from two github repositories
#install_github('rafalcode/GenABEL.data')
#install_github('xamwise/GenABEL')
library('GenABEL')

# Seems to work!!! Let's try an example out

# Fake data
require(GenABEL.data)
data(ge03d2.clean)

gkin <- ibs(ge03d2.clean[,sample(autosomal(ge03d2.clean),1000)], w="freq")
h2ht <- polygenic(height, kin=gkin, ge03d2.clean)

orig <- ge03d2.clean@phdata$height
grres <- h2ht$grresidualY # environmental residuals
pgres <- h2ht$pgresidualY # environmental + predicted breeding values

plot(orig~grres) # changes the axis completely
plot(orig~pgres) # changes the scale completely
plot(grres~pgres) # very similar scale but the values are not identical


