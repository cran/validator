intVal <-
function (y, x, index = "all")
{
    clres <- y
    #### NEW ####
    cluster <- predict(clres, x)
    x <- as.matrix(x)
    clsize <- table(cluster)
    centers <- clres@centers
    #############
     varwithinss <- function(x, centers, cluster) {
        x <- (x - centers[cluster, ])^2
        varwith <- aggregate(x, by=list(cluster), FUN=sum)
        varwithins <- as.matrix(varwith[,-1])
        return(varwithins)
      }
    withinss <- function(varwithins) {
        withins <- apply(varwithins, 1, sum)
        return(withins)
    }
    gss <- function(x, clsize, withins) {
        n <- sum(clsize)
        k <- length(clsize)
        allmean <- apply(x, 2, mean)
        dmean <- sweep(x, 2, allmean, "-")
        allmeandist <- sum(dmean^2)
        wgss <- sum(withins)
        bgss <- allmeandist - wgss
        zgss <- list(wgss = wgss, bgss = bgss)
        return(zgss)
    }
    vargss <- function(x, clsize, varwithins) {
        nvar <- dim(x)[2]
        varallmean <- apply(x, 2, mean)
        vardmean <- (sweep(x, 2, varallmean, "-"))^2
        varallmeandist <- apply(vardmean, 2, sum)
        varwgss <- apply(varwithins, 2, sum)
        vartss <- varallmeandist
        varbgss <- vartss - varwgss
        zvargss <- list(vartss = vartss, varbgss = varbgss)
        return(zvargss)
    }
    ttww <- function(x, clsize, cluster) {
        n <- sum(clsize)
        k <- length(clsize)
        w <- 0
        tt <- cov(x) * (n - 1)
        for (l in 1:k) w <- w + cov(x[cluster == l, ]) * (clsize[l] - 1)
        zttw <- list(tt = tt, w = w)
        return(zttw)
    }
    calinski <- function(zgss, clsize) {
        n <- sum(clsize)
        k <- length(clsize)
        vrc <- (zgss$bgss/(k - 1))/(zgss$wgss/(n - k))
        return(vrc = vrc)
    }
    cindex <- function(withins, minmaxd, clsize) {
        dw <- sum(withins * clsize)
        cindex <- (dw - minmaxd$mindw)/(minmaxd$maxdw - minmaxd$mindw)
        return(cindex)
    }
    db <- function(withins, centers, cluster) {
        mse <- withins/table(cluster)
        r <- outer(mse, mse, "+")/as.matrix(dist(centers, diag = TRUE))
        diag(r) <- 0
        db <- mean(apply(r, 1, max))
        return(db)
    }
    hartigan <- function(zgss) {
        hart <- log(zgss$bgss/zgss$wgss)
        return(hart)
    }
    ratkowsky <- function(zvargss, clsize) {
        k <- length(clsize)
        rat <- mean(sqrt(zvargss$varbgss/zvargss$vartss))
        rat <- rat/sqrt(k)
        return(rat)
    }
    scott <- function(zttw, clsize) {
        n <- sum(clsize)
        dettt <- prod(eigen(zttw$tt)$values)
        detw <- prod(eigen(zttw$w)$values)
        scott <- n * log(dettt/detw)
        return(scott)
    }
    marriot <- function(zttw, clsize) {
        k <- length(clsize)
        detw <- prod(eigen(zttw$w)$values)
        mar <- (k^2) * detw
        return(mar)
    }
    ball <- function(withins, clsize) {
        ball <- sum(withins)/length(clsize)
    }
    tracecovw <- function(zttw) {
        trcovw <- sum(diag(cov(zttw$w)))
        return(trcovw)
    }
    tracew <- function(zttw) {
        tracew <- sum(diag(zttw$w))
        return(tracew)
    }
    friedman <- function(zttw) {
        b <- zttw$tt - zttw$w
        fried <- sum(diag(solve(zttw$w) %*% b))
        return(fried)
    }
    rubin <- function(zttw) {
        dettt <- prod(eigen(zttw$tt)$values)
        detw <- prod(eigen(zttw$w)$values)
        friedm <- dettt/detw
        return(friedm)
    }
    xu <- function(x, clsize, zgss) {
        n <- sum(clsize)
        k <- length(clsize)
        d <- dim(x)[2]
        xuindex <- d * log(sqrt(zgss$wgss/(d * (n^2)))) + log(k)
        return(xuindex)
    }
    varwithins <- varwithinss(x, centers, cluster)
    withins <- withinss(varwithins)
    zgss <- gss(x, clsize, withins)
    zttw <- ttww(x, clsize, cluster)
    index <- pmatch(index, c("calinski", "db", "hartigan",
        "ratkowsky", "scott", "marriot", "ball", "trcovw", "tracew",
        "friedman", "rubin", "xuindex", "dunn", "connectivity", "silhouette",
        "all"))
    if (is.na(index))
        stop("invalid clustering index")
    if (index == -1)
        stop("ambiguous index")
    vecallindex <- numeric(15)
    if (any(index == 1) || (index == 16))
        vecallindex[1] <- calinski(zgss, clsize)
    if (any(index == 2) || (index == 16))
        vecallindex[2] <- db(withins, centers, cluster)
    if (any(index == 3) || (index == 16))
        vecallindex[3] <- hartigan(zgss)
    if (any(index == 4) || (index == 16)) {
        zvargss <- vargss(x, clsize, varwithins)
        vecallindex[4] <- ratkowsky(zvargss, clsize)
      }
    if (any(index == 5) || (index == 16))
        vecallindex[5] <- scott(zttw, clsize)
    if (any(index == 6) || (index == 16))
        vecallindex[6] <- marriot(zttw, clsize)
    if (any(index == 7) || (index == 16))
        vecallindex[7] <- ball(withins, clsize)
    if (any(index == 8) || (index == 16))
        vecallindex[8] <- tracecovw(zttw)
    if (any(index == 9) || (index == 16))
        vecallindex[9] <- tracew(zttw)
    if (any(index == 10) || (index == 16))
        vecallindex[10] <- friedman(zttw)
    if (any(index == 11) || (index == 16))
        vecallindex[11] <- rubin(zttw)
    if (any(index == 12) || (index == 16))
        vecallindex[12] <- xu(x, clsize, zgss)
     if (any(index == 13) || (index == 16)) {
        vecallindex[13] <- dunn(clusters=cluster, Data=x)
      }
     if (any(index == 14) || (index == 16)) {
        vecallindex[14] <- connectivity(clusters=cluster, Data=x)
      }
     if (any(index == 15) || (index == 16)) {
        dist <- dist(x)
        vecallindex[15] <- summary(silhouette(cluster, dist))$si.summary[4]
      }
    names(vecallindex) <- c("calinski", "db", "hartigan",
        "ratkowsky", "scott", "marriot", "ball", "trcovw", "tracew",
        "friedman", "rubin", "xuindex", "dunn", "connectivity", "silhouette")
    if (index < 16)
        vecallindex <- vecallindex[index]
    return(vecallindex)
}

