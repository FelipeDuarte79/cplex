## A linear program (this is lpex1.c in the CPLEX examples)
cvec <- c(1,2,3)
Amat <- matrix(c(-1,1,1,-1,3,-1),byrow=TRUE,nc=3)
bvec <- c(20,-30)
ub <- c(40,Inf,Inf)

res <- Rcplex(cvec,Amat,bvec,ub=ub,objsense="max",sense=c('L','G'))
print(res)

## A linear program with random data
## use the barrier method
n = 20; m = 25
nnz <- trunc(.2 * m * n)

## entries in simple_triplet_matrix clas
##  *must* be in column major order
nnz <- sort(sample(m*n,nnz,replace=FALSE)-1)
Amat <- simple_triplet_matrix(
             i = (nnz %% m) + 1,
             j = trunc(nnz/m) + 1,
             v = rnorm(nnz),
             nrow=m,ncol=n)

x0 <- runif(n)
b <- as.matrix(Amat) %*% x0
cvec <- rnorm(n)

res <- Rcplex(cvec,Amat,b,sense='E',control=list(method=4))
print(res)

## A quadratic problem (this is qpex1.c in the CPLEX examples)
cvec <- c(1,2,3)
Qmat <- matrix(c(-33,6,0,
                  6,-22,11.5,
                  0,11.5,-11),
                byrow=TRUE,
                nc=3)
Amat <- matrix(c(-1,1,1,
                  1,-3,1),
               byrow=TRUE,nc=3)
bvec <- c(20,30)
ub <- c(40,Inf,Inf)

res <- Rcplex(cvec,Amat,bvec,Qmat,ub=ub,objsense="max")
print(res)

## A mixed integer linear program (mipex1.c in the CPLEX examples)
cvec <- c(1,2,3,1)
Amat <- matrix(c(-1,1,1,10,
                  1,-3,1,0,
                  0,1,0,-3.5),
               byrow=TRUE, nc=4)
bvec <- c(20,30,0)
lb <- c(0,0,0,2)
ub <- c(40,Inf,Inf,3)
vtype <- c(rep("C",3),"I")

res <- Rcplex(cvec,Amat,bvec,lb=lb,ub=ub,sense=c("L","L","E"),
              objsense="max",vtype=vtype)
print(res)

## A mixed integer quadratic program
cvec <- c(1,2,3,1)
Qmat <- matrix(c(-33,6,0,0,
                  6,-22,11.5,0,
                  0,11.5,-11,0,
                  0,0,0,0),
               byrow=TRUE, nc=4)
Amat <- matrix(c(-1,1,1,10,
                  1,-3,1,0,
                  0,1,0,-3.5),
               byrow=TRUE, nc=4)
bvec <- c(20,30,0)
ub <- c(40,Inf,Inf,3)
vtype <- c(rep("C",3),"I")

res <- Rcplex(cvec,Amat,bvec,Qmat=Qmat,ub=ub,sense=c("L","L","E"),
              objsense="max",vtype=vtype)
print(res)
Rcplex.close()
