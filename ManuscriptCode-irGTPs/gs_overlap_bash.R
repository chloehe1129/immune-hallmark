
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])

seq = c(seq(1,5200,10)[seq(1,length(seq(1,5200,10)),1)],5200,5219)
x=seq[i]
y=seq[i+1]
source("/rsrch4/home/bcb/she4/hallmark/code/gs-overlap.R")


# for(i in 1:(length(seq)-1)){
#   if(i==521){
#     x=seq[i]
#     y=seq[i+1]
#     source("/rsrch4/home/bcb/she4/hallmark/code/gs-overlap.R")
#   }else{
#     x=seq[i]
#     y=seq[i+1]-1
#     source("/rsrch4/home/bcb/she4/hallmark/code/gs-overlap.R")
#   }
# }

