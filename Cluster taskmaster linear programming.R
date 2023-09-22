# Cluster code optimization
library(Rglpk)

# You need to time simulation speeds for each of the parameter sets and save the timings in a file
timings<-readRDS("D:\\Dropbox\\fileset_144_timings_for_64_final_sets.rds")


N<-256 # (4*64 = 256) tasks because the 64 simulation conditions were each split into 4 computing tasks
M<-72  # Across 72 cores

# See Rglpk for details on how to write the linear programming
item_times<-rep(timings/4,4)

obj<-c(1,rep(0,M + M*N))

mat_1<-Reduce(cbind,list(rep(1,M),diag(-1,nrow=M),matrix(0,ncol=M*N,nrow=M)))

item_time_mat<-matrix(0,nrow=M,ncol=M*N)
for(i in 1:M){item_time_mat[i,seq(N*(i-1)+1,N*(i-1)+length(item_times))]<-item_times}
mat_2<-Reduce(cbind,list(rep(0,M),diag(-1,nrow=M),item_time_mat))

mat_3<-Reduce(cbind,c(list(matrix(0,nrow=N,ncol=1+M)),rep(list(diag(N)),M)))

constraint<-matrix(0,M,M*N)
for(i in 1:M){constraint[i,seq(N*(i-1)+1,N*(i-1)+length(item_times))]<-1}
mat_4<-Reduce(cbind,list(matrix(0,nrow=M,ncol=1+M),constraint))

mat<-Reduce(rbind,list(mat_1,mat_2,mat_3,mat_4))

dir<-c(rep(">=",M),rep("==",M),rep("==",N),rep("==",M))

rhs<-c(rep(0,M),rep(0,M),rep(1,N),sapply(splitIndices(64*4,72),length))

types1<-c(rep("C",1+M),rep("B",M*N))

# Runs the Rglpk solver
solution<-Rglpk_solve_LP(obj, mat, dir, rhs, types=types1, max=FALSE, verbose=TRUE,control=list(tm_limit=7.2e+6)) # Two hour search

# Checking the solution
solution_test<-matrix(solution$solution[74:18505],nrow=M,ncol=N,byrow=T)

rowSums(solution_test)
all(sapply(splitIndices(64*4,72),length)==rowSums(solution_test))

colSums(solution_test)

# Extracting the file sets for the nodes
which(solution_test==1,arr.ind = T) -> assignments
assignments<-assignments[order(assignments[,1]),]

# Checking the result
which(tapply(assignments[,2],assignments[,1],function(x){any(duplicated(x %% 64))}))
cbind(assignments[,1],assignments[,2] %% 64)


assignments<-assignments[,2]
all(assignments[order(assignments)]==1:256)

