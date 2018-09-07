### Sachs data analyzer
load("sachs.rda")
sachs.data <- read.table("sachs.data.txt", header=TRUE)
#sachs.data <- log(sachs.data)
#sachs.data <- sachs.data[1:50,]
sachs.data <- apply(sachs.data, 2, scale, scale=FALSE, center=TRUE) 
sachs.data <- as.data.frame(sachs.data)
sachs.modelstring <- "[PKC|Plcg:PIP2][PKA|PKC][Raf|PKC:PKA][Mek|PKC:PKA:Raf][Erk|Mek:PKA][Akt|Erk:PKA:PIP3][P38|PKC:PKA][Jnk|PKC:PKA][Plcg][PIP3|Plcg][PIP2|Plcg:PIP3]"
g <- model2network(sachs.modelstring)
g <- cpdag(g)

#sachs.pc <- pc.stable(sachs.data, alpha=0.05)
sachs.pc <- rsmax2(sachs.data, restrict='pc.stable', maximize='hc', maximize.args=list(score='bge', iss=3))
sachs.pc.boot <- boot.strength(sachs.data, R = 500, algorithm = 'rsmax2', algorithm.args = list(restrict='pc.stable', maximize='hc', maximize.args=list(score='bge', iss=3))) 
boot.dag <- averaged.network(sachs.pc.boot, threshold=0.5)
boot.cpdag <- cpdag(boot.dag)
rsmax2(sachs.data, restrict='pc.stable', maximize='hc', maximize.args=list(score='bge', iss=3))

sachs.pc.cpdag <- cpdag(sachs.pc)
graphviz.compare(g, sachs.pc.cpdag, main=c("Actual", "PC / Actual Comparison - Skeleton"))
compare(skeleton(g), skeleton(sachs.pc.cpdag))

### Stability selection feature % CDF graph

g <- ggplot(boot, aes(x=boot$strength)) + stat_ecdf(geom="point")
g <- g + scale_y_continuous(name='CDF')
g <- g + scale_x_continuous(name='Stability Selection Strength %')
g <- g + theme(plot.title = element_text(hjust = 0.5))
g <- g + ggtitle("CDF of Stability Selection Feature %")
g

### Sachs distribution of individual elements

sachs.data <- sachs.data[sachs.data$Raf < log(200), ]
sachs.data <- sachs.data[sachs.data$PKC < log(60), ]
sachs.data <- sachs.data[sachs.data$Akt < log(200), ]
g1 <- ggplot(sachs.data, aes(x=sachs.data$Raf)) + geom_histogram(color="black", fill="white", binwidth=0.2) 
g1 <- g1 + scale_x_continuous(name='Expression Level') + scale_y_continuous(name='Count')
g1 <- g1 + ggtitle("RAF") + theme(plot.title = element_text(hjust = 0.5))

g2 <- ggplot(sachs.data, aes(x=sachs.data$PKC)) + geom_histogram(color="black", fill="white", binwidth=0.2) 
g2 <- g2 + scale_x_continuous(name='Expression Level') + scale_y_continuous(name='Count')
g2 <- g2 + ggtitle("PKC") + theme(plot.title = element_text(hjust = 0.5))


g3 <- ggplot(sachs.data, aes(x=sachs.data$Akt)) + geom_histogram(color="black", fill="white", binwidth=0.2) 
g3 <- g3 + scale_x_continuous(name='Expression Level') + scale_y_continuous(name='Count')
g3 <- g3 + ggtitle("AKT") + theme(plot.title = element_text(hjust = 0.5))

grid.arrange(g1, g2, g3, nrow = 1)

### Violin plot of model performance

g <- ggplot(res, aes(x=res$Method, y=res$SHD), log="y") 
g <- g + geom_violin() 
g <- g + scale_y_continuous(name='Structural Hamming Distance', trans='log10')
g <- g + scale_x_discrete(name='Method')
g <- g + geom_boxplot(width=0.1)
g <- g + theme(plot.title = element_text(hjust = 0.5))
g <- g + ggtitle("Comparison of methods for simulated Gaussian DAG inference")
g


### Generate simulated data from networks

library("bnlearn")
load("arth150.rda")
for (i in 1:num_tests) {
  N <- 500
  data <- rbn(bn, n=N)
  fn <- paste('FINAL_LARGE_ARTH_sample_', toString(i), ".csv", sep='')
  write.table(data, fn, append = FALSE, sep = ",", dec = ".",
              row.names = FALSE, col.names = TRUE)
}

### Simulated Gaussian Data

for(i in 1:num_tests) {
  print(i)
  N <- runif(n=1, min=20, max=99)
  D <- runif(n=1, min=50, max=50)
  
  sim <- DAGsim(N, D, 3/(D-1), 0)
  bn.true <- empty.graph(paste(c(1:D)))
  amat(bn.true) <- sim$Adjacency.matrix
  bn.true.cpdag = cpdag(bn.true)
  sim.df <- data.frame(sim$data)
  
  # Run PC algorithm learning
  bn.pc <- pc.stable(sim.df, alpha=0.01)
  nodes(bn.pc) <- nodes(bn.true)
  bn.pc.cpdag <- cpdag(bn.pc)
  results.pc <- shd(bn.true.cpdag, bn.pc.cpdag)
  
  # Run PC algorithm learning
  bn.iamb <- inter.iamb(sim.df, alpha=0.01)
  nodes(bn.iamb) <- nodes(bn.true)
  bn.iamb.cpdag <- cpdag(bn.iamb)
  results.interiamb <- shd(bn.true.cpdag, bn.iamb.cpdag)
  
  # Run HC from NULL graph - 5 restarts, 10 deletes
  bn.hc.null <- hc(sim.df, restart=10, perturb=5, max.iter=100000) 
  nodes(bn.hc.null) <- nodes(bn.true)
  bn.hc.null.cpdag <- cpdag(bn.hc.null)
  results.hc.null <- shd(bn.true.cpdag, bn.hc.null.cpdag)
  
  # Hybrid - HC (same params) + PC skel
  bn.hybrid <- rsmax2(sim.df, restrict='pc.stable', maximize='hc',
                      restrict.args=list(alpha=0.01), 
                      maximize.args=list(restart=10, perturb=5, max.iter=100000))
  nodes(bn.hybrid) <- nodes(bn.true)
  bn.hybrid.cpdag <- cpdag(bn.hybrid)
  results.hybrid <- shd(bn.true.cpdag, bn.hybrid.cpdag)
  
  # Hybrid - HC (same params) + InterIAMB
  bn.hybridiamb <- rsmax2(sim.df, restrict='inter.iamb', maximize='hc',
                          restrict.args=list(alpha=0.01), 
                          maximize.args=list(restart=10, perturb=5, max.iter=100000))
  nodes(bn.hybridiamb) <- nodes(bn.true)
  bn.hybridiamb.cpdag <- cpdag(bn.hybridiamb)
  results.hybridiamb <- shd(bn.true.cpdag, bn.hybridiamb.cpdag)
  
  # HC - Random graph
  hc.init.graph <- random.graph(nodes(bn.true), num=1, prob = 3/(D-1))
  nodes(hc.init.graph) <- colnames(sim.df)
  bn.hc.random <- hc(sim.df, start=hc.init.graph, restart=10, perturb=5, max.iter=100000) 
  nodes(bn.hc.random) <- nodes(bn.true)
  bn.hc.random.cpdag <- cpdag(bn.hc.random)
  results.hc.random <- shd(bn.hc.random.cpdag, bn.true.cpdag)
  
  
  #Empty graph 
  bn.empty <- empty.graph(nodes(bn.true), num=1)
  results.empty <- shd(bn.true.cpdag, bn.empty)
  #graphviz.compare(bn.true.cpdag, bn.hybrid.cpdag)
  
  res[nrow(res)+1, ] = list("PC", results.pc) 
  res[nrow(res)+1, ] = list("HC-Null", results.hc.null)
  res[nrow(res)+1, ] = list("HC-Random", results.hc.random) 
  res[nrow(res)+1, ] = list("Hybrid-PC", results.hybrid)
  res[nrow(res)+1, ] = list("Hybrid-InterIAMB", results.hybridiamb)
  res[nrow(res)+1, ] = list("Empty", results.empty)
  res[nrow(res)+1, ] = list("InterIAMB", results.interiamb)
}

### Network SHD distribuiton calculator

#res <- data.frame(matrix(vector(), 0, 3,
#                       dimnames=list(c(), c("Method", "Network", "SHD"))),
#              stringsAsFactors=T)

load("arth150.rda")
network <- "ARTH-LARGE"
N <- 500
bn.true <- bn
bn.true.cpdag <- cpdag(bn.true)

for(i in 1:num_tests) {
  print(i)
  sim.df <- rbn(bn, n=N)
  # Run PC algorithm learning
  bn.pc <- pc.stable(sim.df, alpha=0.01)
  nodes(bn.pc) <- nodes(bn.true)
  bn.pc.cpdag <- cpdag(bn.pc)
  results.pc <- shd(bn.true.cpdag, bn.pc.cpdag)
  
  # Run InterIAMB algorithm learning
  bn.iamb <- inter.iamb(sim.df, alpha=0.01)
  nodes(bn.iamb) <- nodes(bn.true)
  bn.iamb.cpdag <- cpdag(bn.iamb)
  results.interiamb <- shd(bn.true.cpdag, bn.iamb.cpdag)
  
  # Hybrid - HC (same params) + PC skel
  bn.hybrid <- rsmax2(sim.df, restrict='pc.stable', maximize='hc',
                      restrict.args=list(alpha=0.01), 
                      maximize.args=list(restart=20, perturb=10, max.iter=100000, score='bge', iss=3))
  nodes(bn.hybrid) <- nodes(bn.true)
  bn.hybrid.cpdag <- cpdag(bn.hybrid)
  results.hybrid <- shd(bn.true.cpdag, bn.hybrid.cpdag)
  
  # Hybrid - HC (same params) + InterIAMB
  bn.hybridiamb <- rsmax2(sim.df, restrict='inter.iamb', maximize='hc',
                          restrict.args=list(alpha=0.01), 
                          maximize.args=list(restart=20, perturb=10, max.iter=100000, score='bge', iss=3))
  nodes(bn.hybridiamb) <- nodes(bn.true)
  bn.hybridiamb.cpdag <- cpdag(bn.hybridiamb)
  results.hybridiamb <- shd(bn.true.cpdag, bn.hybridiamb.cpdag)
  
  res[nrow(res)+1, ] = list("PC", network, results.pc) 
  res[nrow(res)+1, ] = list("Hybrid-PC", network, results.hybrid)
  res[nrow(res)+1, ] = list("Hybrid-InterIAMB", network, results.hybridiamb)
  res[nrow(res)+1, ] = list("InterIAMB", network, results.interiamb)
}


### Network SHd distribution plotter

r1 <- res[res$Network == 'ECOLI-SMALL', ]
g1 <- ggplot(r1, aes(x=r1$Method, y=r1$SHD), log="y") 
g1 <- g1 + geom_violin() 
g1 <- g1 + scale_y_continuous(name='Structural Hamming Distance', trans='log10')
g1 <- g1 + scale_x_discrete(name='Method')
g1 <- g1 + geom_boxplot(width=0.1)
g1 <- g1 + theme(plot.title = element_text(hjust = 0.5))
g1 <- g1 + ggtitle("ECOLI70 Inference, N=30")

## ECOLI LARGE

r2 <- res[res$Network == 'ECOLI-LARGE', ]
g2 <- ggplot(r2, aes(x=r2$Method, y=r2$SHD), log="y") 
g2 <- g2 + geom_violin() 
g2 <- g2 + scale_y_continuous(name='Structural Hamming Distance', trans='log10')
g2 <- g2 + scale_x_discrete(name='Method')
g2 <- g2 + geom_boxplot(width=0.1)
g2 <- g2 + theme(plot.title = element_text(hjust = 0.5))
g2 <- g2 + ggtitle("ECOLI70 Inference, N=500")


## NIAB SMALL

r3 <- res[res$Network == 'NIAB-SMALL', ]
g3 <- ggplot(r3, aes(x=r3$Method, y=r3$SHD), log="y") 
g3 <- g3 + geom_violin() 
g3 <- g3 + scale_y_continuous(name='Structural Hamming Distance', trans='log10')
g3 <- g3 + scale_x_discrete(name='Method')
g3 <- g3 + geom_boxplot(width=0.1)
g3 <- g3 + theme(plot.title = element_text(hjust = 0.5))
g3 <- g3 + ggtitle("MAGIC-NIAB Inference, N=100")

## NIABLARGE

r4 <- res[res$Network == 'NIAB-LARGE', ]
g4 <- ggplot(r4, aes(x=r4$Method, y=r4$SHD), log="y") 
g4 <- g4 + geom_violin() 
g4 <- g4 + scale_y_continuous(name='Structural Hamming Distance', trans='log10')
g4 <- g4 + scale_x_discrete(name='Method')
g4 <- g4 + geom_boxplot(width=0.1)
g4 <- g4 + theme(plot.title = element_text(hjust = 0.5))
g4 <- g4 + ggtitle("MAGIC-NIAB Inference, N=500")

## ARTHSMALL

r5 <- res[res$Network == 'ARTH-SMALL', ]
g5 <- ggplot(r5, aes(x=r5$Method, y=r5$SHD), log="y") 
g5 <- g5 + geom_violin()
g5 <- g5 + scale_y_continuous(name='Structural Hamming Distance', trans='log10')
g5 <- g5 + scale_x_discrete(name='Method')
g5 <- g5 + geom_boxplot(width=0.1)
g5 <- g5 + theme(plot.title = element_text(hjust = 0.5))
g5 <- g5 + ggtitle("ARTH150 Inference, N=50")

## ARTH LARGE

r6 <- res[res$Network == 'ARTH-LARGE', ]
g6 <- ggplot(r6, aes(x=r6$Method, y=r6$SHD), log="y") 
g6 <- g6 + geom_violin() 
g6 <- g6 + scale_y_continuous(name='Structural Hamming Distance', trans='log10')
g6 <- g6 + scale_x_discrete(name='Method')
g6 <- g6 + geom_boxplot(width=0.1)
g6 <- g6 + theme(plot.title = element_text(hjust = 0.5))
g6 <- g6 + ggtitle("ARTH150 Inference, N=500")

grid.arrange(g1,g2,g3,g4,g5,g6, nrow = 3)

### Hamming distance recovery for varying networks

#res <- data.frame(matrix(vector(), 0, 4, dimnames=list(c(), c("Method", "Network", "ND", "HD"))), stringsAsFactors=T)

for(t in 1:num_tests){
  print(t)
  load("ecoli70.rda")
  D <- narcs(skeleton(bn))
  Nfracs <- c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
  for (j in 1:length(Nfracs)) {
    N <- as.integer(D*Nfracs[j])
    data <- rbn(bn, n=N)
    
    iamb <- inter.iamb(data, alpha=0.01)
    nodes(iamb) <- nodes(bn)
    hd <- hamming(iamb, skeleton(bn))
    res[nrow(res)+1, ] = list("InterIAMB", "ECOLI", Nfracs[j], hd) 
  }
  load("magic-niab.rda")
  D <- narcs(skeleton(bn))
  Nfracs <- c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
  for (j in 1:length(Nfracs)) {
    N <- as.integer(D*Nfracs[j])
    data <- rbn(bn, n=N)
    
    iamb <- inter.iamb(data, alpha=0.01)
    nodes(iamb) <- nodes(bn)
    hd <- hamming(iamb, skeleton(bn))
    res[nrow(res)+1, ] = list("InterIAMB", "NIAB", Nfracs[j], hd) 
  }
  load("arth150.rda")
  D <- narcs(skeleton(bn))
  Nfracs <- c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
  for (j in 1:length(Nfracs)) {
    N <- as.integer(D*Nfracs[j])
    data <- rbn(bn, n=N)
    
    iamb <- inter.iamb(data, alpha=0.01)
    nodes(iamb) <- nodes(bn)
    hd <- hamming(iamb, skeleton(bn))
    res[nrow(res)+1, ] = list("InterIAMB", "ARTH", Nfracs[j], hd) 
  }
}


### PLOTTER

data.ecoli <- res[res$Network == "ECOLI", ]
data.ecoli <- data.frame(data.ecoli$ND, data.ecoli$HD)
data.ecoli <- aggregate(data.ecoli, by=list(data.ecoli$data.ecoli.ND), FUN = mean)
colnames(data.ecoli) <-  c("ND", "ND", "HD")


data.niab <- res[res$Network == "NIAB", ]
data.niab <- data.frame(data.niab$ND, data.niab$HD)
data.niab <- aggregate(data.niab, by=list(data.niab$data.niab.ND), FUN = mean)
colnames(data.niab) <-  c("ND", "ND", "HD")

data.arth <- res[res$Network == "ARTH", ]
data.arth <- data.frame(data.arth$ND, data.arth$HD)
data.arth <- aggregate(data.arth, by=list(data.arth$data.arth.ND), FUN = mean)
colnames(data.arth) <-  c("ND", "ND", "HD")


plot(data.ecoli$ND, data.ecoli$HD/70, col='red', xlab="Sample Size / Number of Vertices", ylab="Normalised Hamming Distance",
     log='x', ylim=c(0,1.2),
     main="InterIAMB Structure recovery for varying sample size")

lines(data.ecoli$ND, data.ecoli$HD/70, col='red', lwd=2)

points(data.niab$ND, data.niab$HD/66, col='blue')
lines(data.niab$ND, data.niab$HD/66,  col='blue', lwd=2)

points(data.arth$ND, data.arth$HD/150, col='green')
lines(data.arth$ND, data.arth$HD/150,  col='green', lwd=2)

points(data.hyb2$HD, data.hyb2$SHD, col='brown')
lines(data.hyb2$HD, data.hyb2$SHD,  col='brown')

lines(c(-10,100), c(-10,100), lty='dotted')

legend("bottomright", "", c("PC", 'InterIAMB', 'HC-InterIAMB', 'HC-PC'),
       fill=c('red', 'blue', 'green', 'brown' ))

### Compelled edge recovery

#res <- data.frame(matrix(vector(), 0, 6, 
#                         dimnames=list(c(), c("Method", "Network", "ND", "EdgeType",  "Recall", "Precision"))), 
#                  stringsAsFactors=T)
load("arth150.rda")

# Prepare the relevant graphs
bn.cpdag <- cpdag(bn)
### Generate the graph of only the compelled edges 
comp <- empty.graph(nodes(bn.cpdag))
arcs(comp) <- compelled.arcs(bn.cpdag)
undir <- empty.graph(nodes(bn.cpdag))
arcs(undir) <- undirected.arcs(bn.cpdag) 

for(t in 1:num_tests){
  print(t)
  D <- narcs(skeleton(bn))
  Nfracs <- c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
  for (j in 1:length(Nfracs)) {
    N <- as.integer(D*Nfracs[j])
    data <- rbn(bn, n=N)
    pc <- rsmax2(data, restrict='pc.stable', restrict.args=list(alpha=0.01), 
                 maximize='hc', maximize.args=list(score='bge', iss=3, restart=20, perturb=10, max.iter=250000))
    nodes(pc) <- nodes(bn)
    pc.cpdag <- cpdag(pc)
    pc.comp <- empty.graph(nodes(pc.cpdag))
    arcs(pc.comp) <- compelled.arcs(pc.cpdag)
    cf <- compare(comp, pc.comp)
    sens <- cf$tp / narcs(comp)
    prec <- cf$tp / narcs(pc.comp)
    res[nrow(res)+1, ] = list("PC", "ARTH", Nfracs[j], "Compelled", sens, prec) 
    
    # Test Skeleton perf
    cf <- compare(skeleton(bn), skeleton(pc.cpdag))
    sens <- cf$tp / narcs(skeleton(bn))
    prec <- cf$tp / narcs(skeleton(pc.cpdag))
    res[nrow(res)+1, ] = list("PC", "ARTH", Nfracs[j], "Skeleton", sens, prec) 
    
    # Test DAG perf
    cf <- compare(bn.cpdag, pc.cpdag)
    sens <- cf$tp / narcs(bn)
    prec <- cf$tp / narcs(pc.cpdag)
    res[nrow(res)+1, ] = list("PC", "ARTH", Nfracs[j], "DAG", sens, prec) 
    
    # Test Non compelled arcs
    pc.undir <- empty.graph(nodes(pc.cpdag))
    arcs(pc.undir) <- undirected.arcs(pc.cpdag)
    cf <- compare(undir, pc.undir)
    sens <- cf$tp / narcs(undir)
    prec <- cf$tp / narcs(pc.undir)
    res[nrow(res)+1, ] = list("PC", "ARTH", Nfracs[j], "Non-Compelled", sens, prec) 
  }
}


### PLOTTER

data.comp <- res[res$Method == "PC" & res$Network=="ARTH" & res$EdgeType=="Compelled", ]
data.comp <- data.frame(data.comp$ND, data.comp$Recall, data.comp$Precision)
data.comp <- aggregate(data.comp, by=list(data.comp$data.comp.ND), FUN = mean, na.rm=TRUE)
colnames(data.comp) <- c("ND", "ND", "Recall", "Precision")

data.skel <- res[res$Method == "PC" & res$Network=="ARTH" & res$EdgeType=="Skeleton", ]
data.skel <- data.frame(data.skel$ND, data.skel$Recall, data.skel$Precision)
data.skel <- aggregate(data.skel, by=list(data.skel$data.skel.ND), FUN = mean, na.rm=TRUE)
colnames(data.skel) <- c("ND", "ND", "Recall", "Precision")

data.cpdag <- res[res$Method == "PC" & res$Network=="ARTH" & res$EdgeType=="DAG", ]
data.cpdag <- data.frame(data.cpdag$ND, data.cpdag$Recall, data.cpdag$Precision)
data.cpdag <- aggregate(data.cpdag, by=list(data.cpdag$data.cpdag.ND), FUN = mean, na.rm=TRUE)
colnames(data.cpdag) <- c("ND", "ND", "Recall", "Precision")

data.undir <- res[res$Method == "PC" & res$Network=="ARTH" & res$EdgeType=="Non-Compelled", ]
data.undir <- data.frame(data.undir$ND, data.undir$Recall, data.undir$Precision)
data.undir <- aggregate(data.undir, by=list(data.undir$data.undir.ND), FUN = mean, na.rm=TRUE)
colnames(data.undir) <- c("ND", "ND", "Recall", "Precision")


plot(data.comp$Recall, data.comp$Precision, col='red', xlab="Recall",
     ylab="Precision", xlim=c(0,1), ylim=c(0,1), 
     main="Precision - Recall curve for CPDAG inference in ARTH150")

lines(data.comp$Recall, data.comp$Precision, col='red', lwd=2)

points(data.skel$Recall, data.skel$Precision, col='blue')
lines(data.skel$Recall, data.skel$Precision,  col='blue', lwd=2)

points(data.cpdag$Recall, data.cpdag$Precision, col='green')
lines(data.cpdag$Recall, data.cpdag$Precision,  col='green', lwd=2)

points(data.undir$Recall, data.undir$Precision, col='brown')
lines(data.undir$Recall, data.undir$Precision,  col='brown', lwd=2)

lines(c(data.skel[10,3], data.undir[10,3], data.cpdag[10,3], data.comp[10,3]),
      c(data.skel[10,4], data.undir[10,4], data.cpdag[10,4], data.comp[10,4]),
      lty='dotted')
lines(c(data.skel[9,3], data.undir[9,3], data.cpdag[9,3], data.comp[9,3]),
      c(data.skel[9,4], data.undir[9,4], data.cpdag[9,4], data.comp[9,4]),
      lty='dotted')
lines(c(data.skel[8,3], data.undir[8,3], data.cpdag[8,3], data.comp[8,3]),
      c(data.skel[8,4], data.undir[8,4], data.cpdag[8,4], data.comp[8,4]),
      lty='dotted')
lines(c(data.skel[7,3], data.comp[7,3], data.cpdag[7,3], data.undir[7,3]),
      c(data.skel[7,4], data.comp[7,4], data.cpdag[7,4], data.undir[7,4]),
      lty='dotted')
for (j in 1:6){
  lines(c(data.skel[j,3], data.comp[j,3], data.cpdag[j,3], data.undir[j,3]),
        c(data.skel[j,4], data.comp[j,4], data.cpdag[j,4], data.undir[j,4]),
        lty='dotted')
}

#lines(c(-10,100), c(-10,100), lty='dotted')

legend("bottomright", "", c('Skeleton', 'CPDAG-All', "CPDAG-Compelled", 'CPDAG-Undirected'),
       fill=c('blue', 'green', 'red', 'brown' ))

### Network Recovery

#res <- data.frame(matrix(vector(), 0, 4, dimnames=list(c(), c("Method", "Network", "ND", "HD"))), stringsAsFactors=T)

for(t in 1:num_tests){
  print(t)
  load("ecoli70.rda")
  D <- narcs(skeleton(bn))
  Nfracs <- c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
  for (j in 1:length(Nfracs)) {
    N <- as.integer(D*Nfracs[j])
    data <- rbn(bn, n=N)
    
    iamb <- inter.iamb(data, alpha=0.01)
    nodes(iamb) <- nodes(bn)
    hd <- hamming(iamb, skeleton(bn))
    res[nrow(res)+1, ] = list("InterIAMB", "ECOLI", Nfracs[j], hd) 
  }
  load("magic-niab.rda")
  D <- narcs(skeleton(bn))
  Nfracs <- c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
  for (j in 1:length(Nfracs)) {
    N <- as.integer(D*Nfracs[j])
    data <- rbn(bn, n=N)
    
    iamb <- inter.iamb(data, alpha=0.01)
    nodes(iamb) <- nodes(bn)
    hd <- hamming(iamb, skeleton(bn))
    res[nrow(res)+1, ] = list("InterIAMB", "NIAB", Nfracs[j], hd) 
  }
  load("arth150.rda")
  D <- narcs(skeleton(bn))
  Nfracs <- c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
  for (j in 1:length(Nfracs)) {
    N <- as.integer(D*Nfracs[j])
    data <- rbn(bn, n=N)
    
    iamb <- inter.iamb(data, alpha=0.01)
    nodes(iamb) <- nodes(bn)
    hd <- hamming(iamb, skeleton(bn))
    res[nrow(res)+1, ] = list("InterIAMB", "ARTH", Nfracs[j], hd) 
  }
}


### PLOTTER

data.ecoli <- res[res$Network == "ECOLI", ]
data.ecoli <- data.frame(data.ecoli$ND, data.ecoli$HD)
data.ecoli <- aggregate(data.ecoli, by=list(data.ecoli$data.ecoli.ND), FUN = mean)
colnames(data.ecoli) <-  c("ND", "ND", "HD")


data.niab <- res[res$Network == "NIAB", ]
data.niab <- data.frame(data.niab$ND, data.niab$HD)
data.niab <- aggregate(data.niab, by=list(data.niab$data.niab.ND), FUN = mean)
colnames(data.niab) <-  c("ND", "ND", "HD")

data.arth <- res[res$Network == "ARTH", ]
data.arth <- data.frame(data.arth$ND, data.arth$HD)
data.arth <- aggregate(data.arth, by=list(data.arth$data.arth.ND), FUN = mean)
colnames(data.arth) <-  c("ND", "ND", "HD")


plot(data.ecoli$ND, data.ecoli$HD/70, col='red', xlab="Sample Size / Number of Vertices", ylab="Normalised Hamming Distance",
     log='x', ylim=c(0,1.2),
     main="InterIAMB Structure recovery for varying sample size")

lines(data.ecoli$ND, data.ecoli$HD/70, col='red', lwd=2)

points(data.niab$ND, data.niab$HD/66, col='blue')
lines(data.niab$ND, data.niab$HD/66,  col='blue', lwd=2)

points(data.arth$ND, data.arth$HD/150, col='green')
lines(data.arth$ND, data.arth$HD/150,  col='green', lwd=2)


legend("bottomleft", "", c("ECOLI", 'MAGIC-NIAB', 'ARTH'),
       fill=c('red', 'blue', 'green' ))

### Model averaging

#res <- data.frame(matrix(vector(), 0, 6,
#                         dimnames=list(c(), c("Method", "Threshold", "HD", "PredEdges", "Recall", "Prec"))),
#                 stringsAsFactors=T)

load("magic-niab.rda")
for(t in 1:num_tests){
  print(t)
  # Run PC algorithm learning
  N <- as.integer(D*Nfracs[j])
  data <- rbn(bn, n=50)
  pc <- pc.stable(data, alpha=0.05)
  nodes(pc) <- nodes(bn)
  cf <- compare(skeleton(bn), skeleton(pc))
  recall <- cf$tp / (cf$tp + cf$fn)
  precision <- cf$tp / (cf$tp + cf$fp)
  edges <- narcs(pc)
  hd <- hamming(pc, skeleton(bn))
  res[nrow(res)+1, ] = list("PC", 0, hd, edges, recall, precision) 
  
  # Run PC Bootstrapped
  pc.boot <- boot.strength(data, R=200, algorithm = 'pc.stable', algorithm.args = list(alpha=0.05)) 
  thresh <- 0.5
  pc.boot.skel <- skeleton(averaged.network(pc.boot, threshold=thresh))
  cf <- compare(skeleton(bn), pc.boot.skel)
  recall <- cf$tp / (cf$tp + cf$fn)
  precision <- cf$tp / (cf$tp + cf$fp)
  edges <- narcs(pc.boot.skel)
  hd <- hamming(pc.boot.skel, skeleton(bn))
  res[nrow(res)+1, ] = list("PC - Averaged 0.01", thresh, hd, edges, recall, precision) 
}


res <- data.frame(matrix(vector(), 0, 4, dimnames=list(c(), c("Method", "Network", "ND", "HD"))), stringsAsFactors=T)

for(t in 1:num_tests){
  print(t)
  load("ecoli70.rda")
  D <- length(nodes(bn))
  Nfracs <- c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
  for (j in 1:length(Nfracs)) {
    N <- as.integer(D*Nfracs[j])
    data <- rbn(bn, n=N)
    
    iamb <- inter.iamb(data, alpha=0.01)
    nodes(iamb) <- nodes(bn)
    hd <- hamming(iamb, skeleton(bn))
    res[nrow(res)+1, ] = list("InterIAMB", "ECOLI", Nfracs[j], hd) 
  }
  load("magic-niab.rda")
  D <- length(nodes(bn))
  Nfracs <- c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
  for (j in 1:length(Nfracs)) {
    N <- as.integer(D*Nfracs[j])
    data <- rbn(bn, n=N)
    
    iamb <- inter.iamb(data, alpha=0.01)
    nodes(iamb) <- nodes(bn)
    hd <- hamming(iamb, skeleton(bn))
    res[nrow(res)+1, ] = list("InterIAMB", "NIAB", Nfracs[j], hd) 
  }
  load("arth150.rda")
  D <- length(nodes(bn))
  Nfracs <- c(1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
  for (j in 1:length(Nfracs)) {
    N <- as.integer(D*Nfracs[j])
    data <- rbn(bn, n=N)
    
    iamb <- inter.iamb(data, alpha=0.01)
    nodes(iamb) <- nodes(bn)
    hd <- hamming(iamb, skeleton(bn))
    res[nrow(res)+1, ] = list("InterIAMB", "ARTH", Nfracs[j], hd) 
  }
}


### PLOTTER

data.ecoli <- res[res$Network == "ECOLI", ]
data.ecoli <- data.frame(data.ecoli$ND, data.ecoli$HD)
data.ecoli <- aggregate(data.ecoli, by=list(data.ecoli$data.ecoli.ND), FUN = mean)
colnames(data.ecoli) <-  c("ND", "ND", "HD")


data.niab <- res[res$Network == "NIAB", ]
data.niab <- data.frame(data.niab$ND, data.niab$HD)
data.niab <- aggregate(data.niab, by=list(data.niab$data.niab.ND), FUN = mean)
colnames(data.niab) <-  c("ND", "ND", "HD")

data.arth <- res[res$Network == "ARTH", ]
data.arth <- data.frame(data.arth$ND, data.arth$HD)
data.arth <- aggregate(data.arth, by=list(data.arth$data.arth.ND), FUN = mean)
colnames(data.arth) <-  c("ND", "ND", "HD")


plot(data.ecoli$ND, data.ecoli$HD/70, col='red', xlab="Sample Size / Number of Vertices", ylab="Normalised Hamming Distance",
     log='x', ylim=c(0,1.2),
     main="InterIAMB skeleton inference for varying sample size")

lines(data.ecoli$ND, data.ecoli$HD/70, col='red', lwd=2)

points(data.niab$ND, data.niab$HD/66, col='blue')
lines(data.niab$ND, data.niab$HD/66,  col='blue', lwd=2)

points(data.arth$ND, data.arth$HD/150, col='green')
lines(data.arth$ND, data.arth$HD/150,  col='green', lwd=2)


legend("bottomleft", "", c("ECOLI", 'MAGIC-NIAB', 'ARTH'),
       fill=c('red', 'blue', 'green' ))

