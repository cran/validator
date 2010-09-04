extVal <-
function(x, y, index="all") {

x <- as.integer(x)
y <- as.integer(y)

# Creating the similarity table
sim <- std.ext(x, y)
  
# Implementation of the index Hamann (1961) and Hubert (1977)
# Range of index: [-1,1]
Hamann <- function(sim) { 
  ham <- ((sim$SS + sim$DD) - (sim$SD + sim$DS)) / (sim$SS + sim$SD + sim$DS + sim$DD)
  return(ham)
}

# Implementation of the index Czekanowski (1932), Dice (1945), Gower and Legendre (1986)
# Range of index: [0,1]
Czekanowski <- function(sim) {  
  cze <- 2*sim$SS/(2*sim$SS + sim$SD + sim$DS)
  return(cze)
}

# Implementation of the index Kulczynski (1927)
# Range of index: [0,1]
Kulczynski <- function(sim) {
  kul <- 0.5*((sim$SS / (sim$SS + sim$SD)) + (sim$SS / (sim$SS + sim$DS)))
  return(kul)
}

# Implementation of the index McConnaughey (1964)
# Range of index: [-1,1]
McConnaughey <- function(sim) {
  mcconn <- ((sim$SS)^2 + sim$SD*sim$DS) / ((sim$SS + sim$SD) * (sim$SS + sim$DS))
  return(mcconn)
}


# Implementation of the index Peirce (1884)
# Range of index: [-1,1]
Peirce <- function(sim) {
  pei <- (sim$SS*sim$DD - sim$SD*sim$DS) / ((sim$SS + sim$DS) * (sim$SD + sim$DD))
  return(pei)
}

# Implementation of the index Wallace (1) (1983)
# Range of index: [0,1]
Wallace1 <- function(sim) {
  wall1 <- (sim$SS) / (sim$SS + sim$SD)
  return(wall1)
}

# Implementation of the index Wallace (2) (1983)
# Range of index: [0,1]
Wallace2 <- function(sim) {
  wall2 <- (sim$SS) / (sim$SS + sim$DS)
  return(wall2)
}

# Implementation of the index Gamma
# Range of index: [-1,1]
Gamma <- function(sim) {
  gam <- (sim$SS*sim$DD + sim$SD*sim$DS) / sqrt((sim$SS+sim$SD) * (sim$SS+sim$DS) * (sim$DS+sim$DD) * (sim$SD+sim$DD))
  return(gam)
}

# Implementation of the index Sokal and Sneath (1) (1963)
# Range of index: [0,1]
Sokal1 <- function(sim) {
  sok1 <- 0.25 * ((sim$SS / (sim$SS+sim$SD)) + (sim$SS / (sim$SS+sim$DS)) + (sim$DD / (sim$DD+sim$SD)) + (sim$DD / (sim$DD+sim$DS)))
  return(sok1)
}

# Implementation of the index Fager and McGowan (1963)
# Range of index: [-0.5,1]
Fager <- function(sim) {
  fag <- (sim$SS / sqrt((sim$SS+sim$SD) * (sim$SS+sim$DS))) - (0.5 / sqrt(sim$SS+sim$SD))
  return(fag)
}

# Implementation of the index Sokal and Sneath (2) (1963)
# Range of index: [0,1]
Sokal2 <- function(sim) {
  sok2 <- (sim$SS / (sim$SS + 2*(sim$SD+sim$DS)))
  return(sok2)
}

# Implementation of the index Sokal and Sneath (3) (1963), Ochiai (1957)
# Range of index: [0,1]
Sokal3 <- function(sim) {
  sok3 <- ((sim$SS * sim$DD) / sqrt((sim$SS+sim$SD) * (sim$SS+sim$DS) * (sim$DS+sim$DD) * (sim$SD+sim$DD)))
  return(sok3)
}

# Implementation of the index Gowor and Legendre (1986), Sokal and Sneath (1963)
# Range of index: [0,1]
Gower <- function(sim) {
  gow <- ((sim$SS + sim$DD) / (sim$SS + 0.5*(sim$SD + sim$DS) + sim$DD))
  return(gow)
}

# Implementation of the index Roger and Tanimoto (1960)
# Range of index: [0,1]
Roger <- function(sim) {
  rog <- ((sim$SS + sim$DD) / (sim$SS + 2*(sim$SD + sim$DS) + sim$DD))
  return(rog)
}

# Implementation of the index Goodman and Kruskal (1954), Yule (1927)
# Range of index: [0,1]
Kruskal <- function(sim) {
  goo <- ((sim$SS*sim$DD - sim$SD*sim$DS) / (sim$SS*sim$DD + sim$SD*sim$DS))
  return(goo)
}

# Implementation of the index Goodman and Kruskal (1954), Yule (1927)
# Range of index: [0,1]
Pearson <- function(sim) {
  phi <- (sim$SS*sim$DD + sim$SD*sim$DS) / ((sim$SS+sim$SD) * (sim$SS+sim$DS) * (sim$DS+sim$DD) * (sim$SD+sim$DD))
  return(phi)
}

# clv.Rand
#Rand <- clv.Rand(sim)

# clv.Jaccard
#Jaccard <- clv.Jaccard(sim)

# clv.Folkes.Mallows
#Folkes <- clv.Folkes.Mallows(sim)

# clv.Russel.Rao
#Russel <- clv.Russel.Rao(sim)

index <- pmatch(index, c("Hamann", "Czekanowski", "Kulczynski", "McConnaughey", "Peirce", "Wallace1", "Wallace2", "Gamma", "Sokal1", "Fager", "Sokal2", "Sokal3", "Gower", "Roger", "Kruskal", "Pearson", "Rand", "Jaccard", "Folkes", "Russel", "all"))
    if (is.na(index))
        stop("invalid clustering index")
    if (index == -1)
        stop("ambiguous index")
    vecallindex <- numeric(20)
    if (any(index == 1) || (index == 21))
        vecallindex[1] <- Hamann(sim)
    if (any(index == 2) || (index == 21))
        vecallindex[2] <- Czekanowski(sim)
    if (any(index == 3) || (index == 21))
        vecallindex[3] <- Kulczynski(sim)
    if (any(index == 4) || (index == 21))
        vecallindex[4] <- McConnaughey(sim)
    if (any(index == 5) || (index == 21))
        vecallindex[5] <- Peirce(sim)
    if (any(index == 6) || (index == 21))
        vecallindex[6] <- Wallace1(sim)
    if (any(index == 7) || (index == 21))
        vecallindex[7] <- Wallace2(sim)
    if (any(index == 8) || (index == 21))
        vecallindex[8] <- Gamma(sim)
    if (any(index == 9) || (index == 21))
        vecallindex[9] <- Sokal1(sim)
    if (any(index == 10) || (index == 21))
        vecallindex[10] <- Fager(sim)
    if (any(index == 11) || (index == 21))
        vecallindex[11] <- Sokal2(sim)
    if (any(index == 12) || (index == 21))
        vecallindex[12] <- Sokal3(sim)
    if (any(index == 13) || (index == 21))
        vecallindex[13] <- Gower(sim)
    if (any(index == 14) || (index == 21))
        vecallindex[14] <- Roger(sim)
    if (any(index == 15) || (index == 21))
        vecallindex[15] <- Kruskal(sim)
    if (any(index == 16) || (index == 21))
        vecallindex[16] <- Pearson(sim)
    if (any(index == 17) || (index == 21))
        vecallindex[17] <- clv.Rand(sim)
    if (any(index == 18) || (index == 21))
        vecallindex[18] <- clv.Jaccard(sim)
    if (any(index == 19) || (index == 21))
        vecallindex[19] <- clv.Folkes.Mallows(sim)
    if (any(index == 20) || (index == 21))
        vecallindex[20] <- clv.Russel.Rao(sim)
    names(vecallindex) <- c("Hamann", "Czekanowski", "Kulczynski", "McConnaughey", "Peirce", "Wallace1", "Wallace2", "Gamma", "Sokal1", "Fager", "Sokal2", "Sokal3", "Gower", "Roger", "Kruskal", "Pearson", "Rand", "Jaccard", "Folkes", "Russel")
    if (index < 21)
        vecallindex <- vecallindex[index]
    return(vecallindex)

}

