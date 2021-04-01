## START OF SCRIPT ##
## clear work space.
rm(list=ls())

## Read-in the values of Rs(CO2 emission), Ts(soil temperature), and Ms(soil moisture). Make sure to change the address to the 
## one on your computer
## NOTE: YOU MUST MAKE SURE THAT THERE ARE NO MISSING VALUE FOR YOUR VARIABLES. 
input.data  <- read.csv (file="d:/rtest2020/raw_data_2014_s_ontario.csv", header=TRUE, sep = ",")

## create the additional variables sq.Ts and ln.Ma.
input.data$sq.Ts_5cm <- input.data$Ts_5cm * input.data$Ts_5cm 
input.data$ln.Ms_30cm <- log(input.data$Ms_30cm)
## HERE YOU HAVE TO RE-ARRANGE THE COLUMNS OF THE INPUT DATA MATRIX TO BE
## IDENTICAL TO THE SEQUWNCE OF THEIR APPEARANCE IN THE MODEL.
input.data <- input.data[c("Rs","Ts_5cm", "sq.Ts_5cm", "Ms_30cm","ln.Ms_30cm")]

##use the nls() function to estimate the coefficients of the Gaussian-Gamma model,
##using the default algorithm = GaussNewton for the nonlinear least square method.
GG.model.fit <- nls(Rs ~ exp(b0 + b1*Ts_5cm + b2*sq.Ts_5cm + b3*Ms_30cm + b4*ln.Ms_30cm),
                    data=input.data,
                    start = list(b0 = 0, b1 = 0, b2 = 0, b3 = 0, b4 = 0))
## if you want to see the non-robust statistics, then activate the following statement. 
#coef(summary(GG.model.fit))
## NOW YOU HAVE TO PUT THE VALUES OF THE DEPENDENT VARIABLE INTO THE COLUMN VECTOR y.
## NOTE: YOU MUST NOT CHANGE y TO ANYTHING ELSE.
y <- matrix(input.data[,1], ncol=1)
## find the number of column for the matrix x that is to be constructed next.
ncol.x <- ncol(input.data)
## NOW YOU HAVE TO PUT THE VALUES OF THE EXPLANATORY VARIABLES INTO THE MATRIX x, STARING FROM COLUMN 2.
## NOTE: THE FIRST COLUMN OF THE MATRIX x MUST BE FILLED WITH 1.
## NOTE: YOU MUST NOT CHANGE x TO ANYTHING ELSE.
x <- as.matrix(cbind(1,input.data[,2:ncol.x]))
## NOW YOU HAVE TO PUT THE VALUES OF THE ESTIMATED COEFFICIENTS INTO THE COLUMN VECTOR b.col.
## NOTE: YOU MUST NOT CHANGE b.col TO ANYTHING ELSE.
b.col <-  coef(GG.model.fit)

###############################################################################
## Beginning of the Module for generating the robust statistics 
###############################################################################
nobs <-  nrow(y)
# Print the number of observations
nobs
y.mean <-  sum(y) / nobs
## y.hat is a column vector containing the predicted emissions.
y.hat <-  exp(x %*% b.col)
y.hat.mat <- matrix()
der <- matrix(0,nrow=nobs,ncol=ncol.x)
for (j in 1:ncol.x) {
der[,j] <-  x[,j] * y.hat}
dert <- t(der)
# der is the Jacobian
diff.i <-  y - y.hat
diff.0 <-  y - y.mean
ssq <-  sum(diff.i * diff.i)
ssq0  <-  sum(diff.0 * diff.0)
der.rb <- der
for (j in 1:ncol.x) {
der.rb[,j] <-  (der[,j] * diff.i * diff.i)}

dert.rb <- t(der.rb)
xpx.rb <-  dert.rb %*% der
# xpx is the information matrix
xpx <-  dert %*% der
# now we change xpx to the inverse of the information matrix
xpx <-  solve(xpx)
df.ess <-  nobs - ncol(x)
nvar <-  ncol(x) -1
# rsrmsq is the Residual Root Mean Square
rsrmsq <-  sqrt(ssq / df.ess)
# rsrmsq0 is the Residual Root Mean Square of the Null Model
rsrmsq0 <-  sqrt(ssq0 / (nobs -1))
r.square <-  (ssq0 - ssq) / ssq0
# print r-square
r.square
adj.r.square <-  1 - (rsrmsq/rsrmsq0)*(rsrmsq/rsrmsq0)
# print Adjusted R-square
adj.r.square
cv.rb <-  xpx %*% xpx.rb %*% xpx
std.err.rb <-  sqrt(diag(cv.rb))
t.ratio.rb <-  b.col / std.err.rb
p.value.rb <-  2 * (1 - pt(abs(t.ratio.rb), df.ess) )
var.names <- colnames(input.data)
var.names <- replace(var.names, c(1),c("intercept"))
coefficient <- b.col
robust.stat <-  cbind(coefficient, std.err.rb, t.ratio.rb, p.value.rb)
rownames(robust.stat) <- c(var.names)
# Print out the table of robust statistics
robust.stat
################################################################################
## End of the Module for generating the robust statistics 
################################################################################

## keep a csv file of the table of robust statistics
write.csv(robust.stat,file="d:/rtest2020/nonlinear.approach.ontario.csv")
## keep all useful output information in a text file. 
##Note: The sink statements work as a pair. 
sink(file="d:/rtest2020/nonlinear.approach.ontario.txt")
print("GAUSSIAN-GAMMA MODEL FOR TEMPERATE FOREST (SOUTHERN ONTARIO): NONLINEAR APPROACH.")
print(robust.stat)
print(c("sample size=",nobs))
print(c("r.square=",r.square,"adj.r.square=", adj.r.square))
sink()
## no more statement beyond this line.

## END OF R- SCRIPT
