library(matrixcalc)
#data for model 1
H3N2_78<- matrix(data = c(66, 87, 25, 22, 4, 13, 14, 15, 9, 4, 0, 4, 4, 9, 1, 0, 0, 4, 3, 1, 0, 0, 0, 1, 1, 0,0,0,0,0), ncol = 5, byrow=TRUE)
H3N2_81<- matrix(data = c(44, 10, 0, 0, 0, 0, 62, 13, 9, 0, 0, 0, 47, 8, 2, 3, 0, 0, 38, 11, 7, 5, 1, 0, 9, 5, 3, 1, 0, 1), ncol = 5)

H3N2_78_1 <- sweep(H3N2_78,2,colSums(H3N2_78),`/`)
H3N2_81_1 <- sweep(H3N2_81,2,colSums(H3N2_81),`/`)

#data for model 2
B_76 <- matrix(data= c(9, 1, 0, 0, 0, 0, 12, 6, 2, 0, 0, 0, 18, 6, 3, 1, 0, 0, 9, 4, 4, 3, 0, 0, 4, 3, 0, 2, 0, 0), ncol=5)
H1N1_79 <- matrix(data= c(15, 11, 0, 0, 0, 0, 12, 17, 21, 0, 0, 0, 4, 4, 4, 5, 0, 0),ncol = 3)

B_76_1 <- sweep(B_76,2,colSums(B_76),`/`)
H1N1_79_1 <- sweep(H1N1_79,2,colSums(H1N1_79),`/`)

set.seed(1)

Data_generator_helper <- function(qlist, susc, index){
  w <- matrix(0, nrow = 6, ncol = susc)
  w[1,] <- sapply(1:susc, function(n,theta){theta^n}, qlist[index]) #populating 1st row w0s: w01,w02,w03,w04,w05
  w[2,1] <- 1 - w_1[1,1] #populating w11
  
  for(s in 2:susc){ #populating rows
     for(j in 2:6){
       if(j<=s){
         w[j,s] <- choose(s,j-1)*w[j,j-1]*(qlist[2] * (qlist[1]^(j-1)))^(s-j+1)
       }
       else{
         w[j,s] <- 1-sum(w[,s])
         break
       }
     }
  }
  return(w)
}

Data_generator <- function(qlist, susc1, susc2) {
  
  w1 <- Data_generator_helper(qlist, susc1, 2)
  w2 <- Data_generator_helper(qlist, susc2, 4)

  return(list("w1" = w1, "w2" = w2))
}

frobenious <- function(val){
  return(sqrt(sum(abs(val)^2, na.rm = TRUE)))
}

Distance <- function(data1, data2, generated_data) {
  #print(generated_data$w1)
  return((1/2) * (frobenious(data1 - generated_data$w1) + frobenious(data2 - generated_data$w2)))
}

ABC <- function(epsilon, n_samples, data1, data2) {
  parameters <- list()
  i <- 1
  for(j in 1 : n_samples) {
    qlist <- runif(4) #qc1 qh1 qc2 qh2
    data <- Data_generator(qlist,ncol(data1),ncol(data2))
    distance <- Distance(data1, data2, data)
    #print(distance)
    if(distance <= epsilon){
      parameters[[i]] <- qlist
      i = i +1
    }
  }
  return(parameters)
}


m1 <- ABC(0.3, 1000, H3N2_78_1, H3N2_81_1)

m2 <- ABC(0.3, 1000, B_76_1, H1N1_79_1)

m1_1 <- as.data.frame(m1)
m1_1_1 <- t(m1_1)
rownames(m1_1_1) <- NULL
colnames(m1_1_1) <- c("qc1","qh1","qc2","qh2")

m2_1 <- as.data.frame(m2)
m2_1_1 <- t(m2_1)
rownames(m2_1_1) <- NULL
colnames(m2_1_1) <- c("qc1","qh1","qc2","qh2")

library(ggplot2)
plot(m1_1_1[,2], m1_1_1[,1], col=c("red"))
plot(m1_1_1[,4], m1_1_1[,3], col=c("blue"))

plot(m2_1_1[,2], m2_1_1[,1], col=c("red"))
plot(m2_1_1[,4], m2_1_1[,3], col=c("blue"))

ggplot(df=data.frame(m1_1_1), aes(x = qh1, y=qc1)) +geom_point(color ="blue") + geom_point(df = data.frame(m2_1_1), aes(x=qh2, y=qc2), color = "red")

ggplot(data=m2_1_1) +
  geom_point(mapping=aes(x=qh1, y=qc1)) +
  geom_point(mapping=aes(x=qh2, y=qc2))
