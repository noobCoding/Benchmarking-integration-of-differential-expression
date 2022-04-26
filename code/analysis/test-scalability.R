for (st in c("nCell", "nGene", "Size", "sqrtSize")) {
  N <- as.numeric(unlist(dfmeta[st]))
  pN <- c()
  for (method in df$X){
    # print(method)
    time <- as.numeric(df[df$X==method,][2:4])
    lrTime = lm(time ~ N)
    pN <- c(pN, lrTime$coefficients[2])
    summary(lrTime)
  }
  df[st] <- pN
}
