fp.from.bstring <- function(fp) {
    size <- length(fp)
    fp <- strsplit(fp,'')[[1]]

    one <- which(fp == '1')
    zero <- which(fp == '0')
    if (length(one) + length(zero) != length(fp)) {
        stop("String can only contain 1 or 0")
    }
    which(fp == '1')
}

fp.and <- function(fp1, fp2, size=1024) {
    smaller <- NULL
    bigger <- NULL
    if (length(fp1) > length(fp2)) {
        smaller <- fp2
        bigger <- fp1
    } else {
        smaller <- fp1
        bigger <- fp2
    }
    i <- bigger %in% smaller
    bigger[i]
}

fp.or <- function(fp1, fp2, size=1024) {
    # a FP structure which is fp1 OR fp2
    r  <- unique(c(fp1,fp2))
    r[order(r)]
}

fp.not <- function(fp1, size=1024) {
    r <- 1:size
    r[ -fp1 ]
}

fp.xor <- function(fp1, fp2, size=1024) {
    tmp1 <- rep(FALSE, size)
    tmp2 <- rep(FALSE, size)
    tmp1[fp1] <- TRUE
    tmp2[fp2] <- TRUE
    tmp3 <- xor(tmp1,tmp2)
    which(tmp3)
}

fp.num <- function(fp) {
    length(fp)
}

#################################

fpNOr <- function( fp1, fp2, size = 1024 ) {
    # num bits on in fp1 but not on in fp2
    n <- length(fp1) - sum(fp1 %in% fp2)
    n
}

fpNAnd <- function( fp1, fp2, size = 1024 ) {
    # num bits on in fp1 and fp2
    n <- sum(fp1 %in% fp2)
    n
}

fpNOff <- function( fp1, fp2, size=1024 ) {
    # num bits off in both fp1 and fp2

    b <- 1:size
    b1 <- b[-fp1] # bit positions off
    b2 <- b[-fp2]
    fpNAnd(b1,b2)
}

# fpRead reads in fingerprints from the format output
# by the CDK fingerprint generator. The return value is
# a list in which each element is a vector of the positions
# in the bit string that are 1.
# Thus the return value does not indicate the length of the 
# original bitstring.

cdk.lf <- function(line) {
    p <- regexpr("{([0-9,\\s]*)}",line,perl=T)
    s <- gsub(',','',substr(line, p+1, p+attr(p,"match.length")-2))
    s <- lapply( strsplit(s,' '), as.numeric )
    s[[1]]
}

moe.lf <- function(line) {
    p <- regexpr("\"([0-9\\s]*)\"",line, perl=T)
    s <- substr(line, p+1, p+attr(p,"match.length")-2)
    s <- lapply( strsplit(s,' '), as.numeric )
    s[[1]]
}

fp.read <- function(f='fingerprint.txt', lf=cdk.lf, header=FALSE) {
    fplist <- list()
    fcon <- file(description=f,open='r')
    lines = readLines(fcon,n=-1)
    if (header) lines = lines[-1]
    c = 1
    for (line in lines) {
        fplist[[c]] <- lf(line)
        c <- c+1
    }
    close(fcon)
    fplist
}

# Need to supply the length of the bit string since fp.read does
# not provide that information
fp.read.to.matrix <- function(f='fingerprint.txt', size=1024, lf=cdk.lf, header=FALSE) {
    fplist <- fp.read(f, lf, header)
    fpmat <- fp.to.matrix(fplist,size)
    fpmat
}
    
fp.distance <- function(fp1,fp2, size=1024, type='tanimoto', ...) {
    c <- fpNAnd(fp1,fp2)
    a <- fpNOr(fp1,fp2)
    b <- fpNOr(fp2,fp1)
    d <- fpNOff(fp1,fp2)

    dist <- NULL

    if (type == 'tanimoto') {
        dist <- c / (a+b+c)
    } else if (type == 'euclidean') {
        dist <- sqrt((d+c) / (a+b+c+d))
    } else if (type == 'dice') {
        dist <- c / (.5*a + .5*b + c)
    } else if (type == 'mt') {
        t1 <- c/(size-d)
        t0 <- d/(size-c)
        phat <- ((size-d) + c)/(2*size)
        dist <- (2-phat)*t1/3 + (1+phat)*t0/3
    }

    dist
}

.readTSV <- function(f="tsv.txt", nmol=1) {

    fcon <- file(description=f,open='r')
    sim <- matrix(1, nr=nmol, nc=nmol)
    for (i in 1:(nmol*(nmol-1)/2)) {
        line <- readLines(fcon,n=1)
        line <- strsplit(line,split=' ')
        m1 <- as.numeric(line[[1]][1]) + 1
        m2 <- as.numeric(line[[1]][2]) + 1
        val <- as.numeric(line[[1]][3])
        sim[m1,m2] <- val
        sim[m2,m1] <- val
    }
    close(fcon)
    sim
}

fp.sim.matrix <- function(fplist, size=1024, type='tanimoto') {
    sim <- matrix(0,nr=length(fplist), nc=length(fplist))
    for (i in 1:(length(fplist)-1)) {
        v <- unlist(lapply( fplist[(i+1):length(fplist)], fp.distance, fp2=fplist[[i]], size=size, type=type))
        sim[i,(i+1):length(fplist)] <- v
        sim[(i+1):length(fplist),i] <- v
    }
    diag(sim) <- 1.0
    sim
}

# Takes the fingerprints, P bits,  for a set of N molecules supplied as
# a list structure and creates an N x P matrix
fp.to.matrix <- function( fplist, size=1024 ) {
    m <- matrix(0, nr=length(fplist), nc=size)
    cnt <- 1
    for ( i in fplist ) {
        m[cnt,i] <- 1
        cnt <- cnt + 1
    }
    m
}

fp.factor.matrix <- function( fplist, size=1024 ) {
    m <- data.frame(fp.to.matrix(fplist,size))
    m[] <- lapply(m, factor, levels=0:1)
    m
}

fp.to.string <- function(fp, size=1024) {
    s <- numeric(size)
    s[fp] <- 1
    paste(s,sep='',collapse='')
}

fp.to.vector <- function(fp, size=1024) {
    coord <- numeric(size)
    coord[fp] <- 1.0 / sqrt(length(fp))
    coord    
}

fp.fold <- function(fp, size=1024) {
    if (size %% 2 != 0) {
        stop('Need to supply a fingerprint of even numbered length')
    }
    bfp <- rep(FALSE, size)
    bfp[fp] <- TRUE
    subfplen <- size/2

    subfp1 <- bfp[1:subfplen]
    subfp2 <- bfp[(subfplen+1):size]

    foldedfp <- fp.xor(subfp1,subfp2, size=subfplen)
    foldedfp
}
