setwd("~/Documents/MSDS16/ImageProcessing/Experiment2")

points <- read.csv('points.csv', header = F)
names(points) <- c('x','y','corr','template','frame','u')
p <- points[,c(1,2,4)]
p <- as.data.frame(scale(p))
#Calculate the dist matrix
d <- dist(p, method = "euclidean")

# Heirarchial 5 groups
p_hh <- hclust(d, method = "ward.D2")
plot(p_hh)
kc = round(sqrt(nrow(points)/2))
phh_sub <- cutree(p_hh, k=kc)
points$phh_sub <- phh_sub
write.csv(points, 'finalPoints.csv')

(points[points$phh_sub > 2,c('x','y', 'phh_sub')])
table(points$phh_sub)


library(fpc)
maxk <- 20  # arbitrary here, you can set this to whatever you like
estimatedK <- pamk(dist(DATA), krange=1:maxk)$nc