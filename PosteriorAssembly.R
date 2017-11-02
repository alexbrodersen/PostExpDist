################################################################################
###
### Code to set up the posterior dstribution
###
################################################################################

rm(list=ls())
library(lpSolve)
library(Rcpp)
setwd("~/Documents/Projects/PostExpDist/")

sourceCpp("FisherInfo.cpp")
sourceCpp("Probs3plm.cpp")
sourceCpp("MLScoring.cpp")

nItems <- 2000
nForms <- 5
nMidpoints <- 20
mLength <- 40

mLength/nItems

1/(nForms*nMidpoints)
kl <- 1
ku <- 2
kl/(nForms*nMidpoints)
ku/(nForms*nMidpoints)

a <- rlnorm(nItems,0,.2)
b <- rnorm(nItems)
c <- rep(0,nItems)

edges <- qnorm(seq(0,1,length.out = (nMidpoints + 1))[2:(nMidpoints)])

lb <- -4
ub <- 4

CalcMidpoints <- function(points,lb,ub){
    x <- c(lb,points,ub)
    y <- NULL
    for(i in 1:(length(x) - 1)){
      y[i] <- (x[i] + x[(i + 1)])/2        
    }
    return(y)
}

Midpoints <- CalcMidpoints(edges,lb,ub)
TRUE
x <- NULL
for(i in 1:nMidpoints){
x <- c(x,rep(FisherInfo(theta = midPoints[i],a = a, b = b, c = c),nForms))
}



makeLengthCon <- function(nItems,nForms,nMidpoints){
    mat <- matrix(0,nrow = nForms*nMidpoints, ncol = nItems*nForms*nMidpoints)
    for(i in 1:(nForms*nMidpoints)){
        mat[i,(1 + (i - 1)*nItems):(i*nItems)] <- 1
    }
    return(mat)
}

makeItemLimits <- function(nItems,nForms,nMidpoints){
    mat <- NULL
    for(i in 1:(nForms*nMidpoints)){
        mat <- cbind(mat,diag(nItems))
    }
    return(mat)
}

f.obj <- x
out <- makeItemLimits(nItems,nForms,nMidpoints)
out2 <- makeLengthCon(nItems,nForms,nMidpoints)


f.con <- rbind(out2,out,out)
f.rhs <- c(rep(mLength,nForms*nMidpoints),rep(3,nItems),rep(1,nItems))
f.dir <- c(rep("=",nForms*nMidpoints),rep("<=",nItems),rep(">=",nItems))
dim(f.con)
length(f.obj)
length(f.dir)
length(f.rhs)
result <- lp("max",f.obj,f.con,f.dir,f.rhs,all.bin = TRUE)
result

rrr <- trimSolution(result,nItems,nForms,nMidpoints)
table(unlist(rrr))



trimSolution <- function(res,nItems,nForms,nMidpoints,trace = F){
    testList <- list()

    for(i in 1:nMidpoints){
        testList[[i]] <- list()
    }

    for(i in 1:nMidpoints){
        for(j in 1:nForms){
            if(trace == T) cat(counter / nForms*nMidpoints)
            counter <- j + nForms*(i - 1)
            testList[[i]][[j]] <- which(res$solution[(1 + (counter - 1)*nItems):(counter*nItems)] == 1)
        }
    }
    return(testList)
}


################################################################################
###
### Big Shadow Test Method
###
################################################################################

x <- NULL

for(i in 1:nMidpoints){
x <- c(x,rep(FisherInfo(theta = midPoints[i],a = a, b = b, c = c),1))
}

length(x)

makeLengthConOne <- function(nItems,nForms,nMidpoints){
    mat <- matrix(0,nrow = nMidpoints, ncol = nItems*nMidpoints)
    for(i in 1:(nMidpoints)){
        mat[i,(1 + (i - 1)*nItems):(i*nItems)] <- 1
    }
    return(mat)
}


makeItemLimitsOne <- function(nItems,nForms,nMidpoints){
    mat <- NULL
    for(i in 1:(nMidpoints)){
        mat <- cbind(mat,diag(nItems))
    }
    return(mat)
}

makeLengthConTwo <- function(nItems,nForms){
    mat <- matrix(0,nrow = nForms, ncol = nItems*nForms)
    for(i in 1:(nForms)){
        mat[i,(1 + (i - 1)*nItems):(i*nItems)] <- 1
    }
    return(mat)
}

makeItemLimitsTwo <- function(nItems,nForms){
    mat <- NULL
    for(i in 1:(nForms)){
        mat <- cbind(mat,diag(nItems))
    }
    return(mat)
}

DivideRes <- function(res,nMidpoints,nItems,nForms,mLength){
    testList <- NULL
    for(i in 1:nMidpoints){
        testList[[i]] <- list()
    }
    
    for(i in 1:nMidpoints){
        cat(i/nMidpoints*100,"%\n")
        f.rhs <- c(rep(mLength,nForms),
                   res$solution[(1 + (i-1)*nItems):(nItems + (i-1)*nItems)])
        f.obj <- rep(res$objective[(1 + (i-1)*nItems):(nItems + (i-1)*nItems)],nForms)
        out <- makeItemLimitsTwo(nItems,nForms)
        out2 <- makeLengthConTwo(nItems,nForms)
        f.con <- rbind(out2,out)
        f.dir <- c(rep("==",nForms),rep("==",nItems))
        dim(f.con)
        length(f.obj)
        length(f.dir)
        length(f.rhs)
        sols <- lp("max",f.obj,f.con,f.dir,f.rhs,all.bin = T)
        testList[[i]] <- trimSolution(sols,nItems,nForms,1)[[1]]
    }
    return(testList)
        
}

gc()
rm(f.con)
rm(out)
rm(out2)
out <- makeItemLimitsOne(nItems,nForms,nMidpoints)
out2 <- makeLengthConOne(nItems,nForms,nMidpoints)

f.obj <- x
f.con <- rbind(out2,out,out)
f.rhs <- c(rep(mLength*nForms,nMidpoints),rep(ku,nItems),rep(kl,nItems))
f.dir <- c(rep("=",nMidpoints),rep("<=",nItems),rep(">=",nItems))
dim(f.con)
length(f.obj)
length(f.dir)
length(f.rhs)
res <- lp("max",f.obj,f.con,f.dir,f.rhs,all.int = TRUE)
res
length(res$solution)/nItems
result$solution
i <- 1
res$constraints
out <- DivideRes(res,nMidpoints,nItems,nForms,mLength)
rrr <- trimSolution(result,nItems,nForms,nMidpoints)
unlist(lapply(rrr,function(x) lapply(x,function(z){
    sum(FisherInfo(0,a[z],b[z],c[z]))
    })))
sort(unlist(out[[1]])) == sort(unlist(rrr[[1]]))
unlist(lapply(out,function(x) lapply(x,function(z){
    sum(FisherInfo(0,a[z],b[z],c[z]))
})))

################################################################################
###
### Bigger Shadow Test Method
###
################################################################################

## The basic idea with the bigger shadow test method is to first partition
## the ability range into a set of windows that are combinations of the smaller
## individual windows. Then, we solve the optimization problems for the larger
## windows. Then, we solve the optimization problem for each (smaller) window.
## Finally, we assembly the tests for each of the smaller problems.
divFactor <- 5



makeLengthConLOne <- function(nItems,nForms,nMidpoints,divFactor){
    if(nMidpoints %% divFactor != 0){
        warning("The modulus of windows and factors is not zero.")
        return(NULL)
    }
    mat <- matrix(0,nrow = nMidpoints/divFactor, ncol = nItems*(nMidpoints/divFactor))
    for(i in 1:(nMidpoints/divFactor)){
        mat[i,(1 + (i - 1)*nItems):(i*nItems)] <- 1
    }
    return(mat)
}

makeItemLimitsLOne <- function(nItems,nForms,nMidpoints,divFactor){
    if(nMidpoints %% divFactor != 0){
        warning("The modulus of windows and factors is not zero.")
        return(NULL)
    }
    mat <- NULL
    for(i in 1:(nMidpoints/divFactor)){
        mat <- cbind(mat,diag(nItems))
        
    }
    return(mat)
}

out <- makeLengthConLOne(nItems,nForms,nMidpoints,divFactor)
out2 <- makeItemLimitsLOne(nItems,nForms,nMidpoints,divFactor)

f.obj <- makeBRSTfOBJ(nItems,nForms,nMidpoints,divFactor,midPoints,a,b,c)
dim(out2)
dim(out)
nMidpoints/divFactor
f.con <- rbind(out,out2,out2)
dim(f.con)
kl <- 1
ku <- 2
mLength*nForms*(nMidpoints/divFactor)
head(f.con)
f.rhs <- c(rep(mLength*nForms*(nMidpoints/divFactor),nMidpoints/divFactor),rep(ku,nItems),rep(kl,nItems))

f.dir <- c(rep("=",nMidpoints/divFactor),rep("<=",nItems),rep(">=",nItems))
length(f.dir)
length(f.rhs)
dim(f.con)
length(f.obj)

out <- lp("max",f.obj,f.con,f.dir,f.rhs,all.int = T)
out
sum(out$solution)
apply(matrix(out$solution,nrow = 4,byrow = T),2,sum)

makeBRSTfOBJ <- function(nItems,nForms,nMidpoints,divFactor,midPoints,a,b,c){
    x <- NULL
    nBiggerWindows <- nMidpoints/divFactor
    counter <- 0
    for(i in 1:nBiggerWindows){
        y <- rep(0,nItems)
        for(j in 1:divFactor){
            counter <- counter + 1
            y <- y + FisherInfo(theta = midPoints[counter],a = a, b = b, c = c)
        }
        x <- c(x,y)
    }
    return(x)
}


makeLengthConLTwo <- function(nItems,nForms){
    mat <- matrix(0,nrow = nForms, ncol = nItems*nForms)
    for(i in 1:(nForms)){
        mat[i,(1 + (i - 1)*nItems):(i*nItems)] <- 1
    }
    return(mat)
}

makeItemLimitsLTwo <- function(nItems,nForms){
    mat <- NULL
    for(i in 1:(nForms)){
        mat <- cbind(mat,diag(nItems))
    }
    return(mat)
}


################################################################################
###
### Run CAT
###
################################################################################

MAP <- function(irvs,a,b,max.it=50,HPM =  0, HPVAR = 1,initTheta = 0){
    theta <- initTheta
    U <- irvs # response vector
    n.items <- length(U)    
    BIGT <- .5 # delta
    for(i in 1:max.it){
        sum.dem <- 0
        sum.num <- 0
        for(j in 1:n.items){
            deviation <- a[j]*(theta-b[j])
            p.hat <- 1/(1 + exp(-deviation))
            WIJ <- p.hat*(1-p.hat)
            VIJ <- (U[j]-p.hat)
            sum.num <- sum.num + a[j]*VIJ
            sum.dem <- sum.dem + a[j]*a[j]*WIJ
        }
        sum.num <- sum.num - (theta-HPM)/HPVAR
        sum.dem <- -sum.dem-(1/HPVAR)
        delta <- sum.num/sum.dem
        if(abs(delta)>BIGT){
            if(delta > 0){
                delta <- BIGT
            } else {
                delta <- -BIGT
            }
        }
        theta <- theta - delta
        if(abs(delta) < .01){
            break
        }
    }
    return(theta)
}


CatSim <- function(thetas,a,b,c,mLength,edges,nForms,rrr,method,trace = T,mod = 100){
    n <- length(thetas)
    theta.hat <- rep(0,n)
    edges <- c(-Inf,edges,Inf)
    I <- matrix(NA,n,mLength)
    U <- matrix(NA,n,mLength)
    theta.out <- matrix(NA,n,mLength)
    for(i in 1:n){
        if(trace == T){
            if(i %% mod == 0){
                cat(i/n*100,"%\n")
            }
        }
        for(j in 1:mLength){

            if(method == "MaxInfo"){                
                info <- FisherInfo(theta.hat[i],a,b,c)
                infoMat <- cbind(t(info),1:length(info))
                items <- na.omit(I[i,])
                if(length(items) != 0){
                    infoMat <- infoMat[-items,]
                }
                infoMat <- infoMat[order(infoMat[,1],decreasing = T),]
                selectedItem <- infoMat[1,2]
                p <- Probs3plm(thetas = thetas[i],a[selectedItem],
                               b[selectedItem],
                               c[selectedItem])
                I[i,j] <- selectedItem
                U[i,j] <- sample(0:1,1,prob = c((1-p),p))
                irvs <- na.omit(U[i,])
                items <- na.omit(I[i,])
                theta.temp <- MAP(irvs,a[items],b[items])
                theta.out[i,j] <- theta.temp
                theta.hat[i] <- theta.temp
            }
            if(method == "PostPool"){
                if(j %in% 1:5){
                selectedItem <- sample(nItems,1)
                p <- Probs3plm(thetas = thetas[i],a[selectedItem],
                               b[selectedItem],
                               c[selectedItem])
                I[i,j] <- selectedItem
                U[i,j] <- sample(0:1,1,prob = c((1-p),p))
                irvs <- na.omit(U[i,])
                items <- na.omit(I[i,])
                theta.temp <- MAP(irvs,a[items],b[items])
                theta.out[i,j] <- theta.temp
                theta.hat[i] <- theta.temp
                next
                }
                whichWindow <- findInterval(theta.hat[i],edges)
                whichForm <- sample(nForms,1)
                windowItems <- rrr[[whichWindow]][[whichForm]]
                availableItems <- windowItems[which(!windowItems %in% I[i,j])]
                selectedItem <- sample(availableItems,1)
                p <- Probs3plm(thetas = thetas[i],a[selectedItem],
                               b[selectedItem],
                               c[selectedItem])
                I[i,j] <- selectedItem
                U[i,j] <- sample(0:1,1,prob = c((1-p),p))
                irvs <- na.omit(U[i,])
                items <- na.omit(I[i,])
                theta.temp <- MAP(irvs,a[items],b[items])
                theta.out[i,j] <- theta.temp
                theta.hat[i] <- theta.temp
            }
        }        
    }
    return(list(theta.out = theta.out,I = I, U = U))
}



n <- 50000
val <- -1
thetas <- rep(val,n)
thetas <- rnorm(n)
res1 <- CatSim(thetas,a,b,c,mLength,edges,nForms,rrr,method = "MaxInfo")
res2 <- CatSim(thetas,a,b,c,mLength,edges,nForms,out,method = "PostPool")
png("ItemExposure.png")
barplot(table(as.vector(res1$I))/n,ylim = c(0,1),main = "Item Exposure Frequency")
dev.off()

mean(res1$theta.out[,mLength] - thetas)
length(names(table(res1$I)))
max(table(res1$I)/n)
hist(res1$theta.out[,mLength])
sqrt(mean((res1$theta.out[,mLength] - thetas)^2))

png("ItemExposure.png")
barplot(table(as.vector(res2$I))/n,ylim = c(0,1),main = "Item Exposure Frequency")
dev.off()

mean(res2$theta.out[,mLength] - thetas)
length(names(table(res2$I)))
hist(table(res2$I)/n,breaks = 50)
max(table(res2$I)/n)
min(table(res2$I)/n)
var(table(res2$I)/n)
quantile(table(res2$I)/n)
hist(res2$theta.out[,mLength] - thetas)
sqrt(mean((res2$theta.out[,mLength] - thetas)^2))

rnorm(10)
for(i in 1:5000){
    if(i %% 100 == 0){
        Sys.sleep(.5)
        plot(res2$theta.out[i,],type = "l",ylim = c(-3,3))
        abline(h = thetas[i])
    }
}


################################################################################
###
### Gurobi Code
###
################################################################################

library(gurobi)
model <- list()
?gurobi
f.rhs <- c(rep(mLength,nForms*nMidpoints),rep(,nItems),rep(25,nItems))
model$A          <- f.con
model$obj        <- f.obj
model$modelsense <- "max"
model$rhs        <- f.rhs
model$sense      <- f.dir
model$vtype      <- 'B'

params <- list(OutputFlag=0)
?gurobi
result <- gurobi(model, params)

print('Solution:')
result$status
print(result$objval)
print(result$x)
