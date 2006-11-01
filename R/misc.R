setGeneric("fold", function(fp) standardGeneric("fold"))
setMethod("fold", "fingerprint",
          function(fp) {
            size <- fp@nbit
            if (size %% 2 != 0) {
              stop('Need to supply a fingerprint of even numbered length')
            }
            bfp <- rep(FALSE, size)
            bfp[fp@bits] <- TRUE

            subfplen <- size/2
            
            b1 <- which(bfp[1:subfplen])
            b2 <- which(bfp[(subfplen+1):size])
            
            subfp1 <- new("fingerprint",
                          nbit=subfplen,
                          bits=b1,
                          provider="R");
            
            subfp2 <- new("fingerprint",
                          nbit=subfplen,
                          bits=b2,
                          provider="R")
            foldedfp <- xor(subfp1,subfp2)
            foldedfp@folded <- TRUE
            return(foldedfp)
          })

setGeneric("euc.vector", function(fp) standardGeneric("euc.vector"))
setMethod("euc.vector", "fingerprint",
          function(fp) {
            coord <- rep(0,length(fp))
            coord[fp@bits] <- 1.0 / sqrt(length(fp))
            coord
          })


setGeneric("distance", function(fp1,fp2,method) standardGeneric("distance"))
setMethod("distance", c("fingerprint", "fingerprint", "missing"),
          function(fp1,fp2) {
            distance(fp1,fp2,"tanimoto")
          })
setMethod("distance", c("fingerprint", "fingerprint", "character"),
          function(fp1,fp2, method=c('tanimoto', 'euclidean', 'mt')) {
            if ( length(fp1) != length(fp2))
              stop("Fingerprints must of the same bit length")
            
            method <- match.arg(method)
            size <- length(fp1)
            
            tmp <- fp1 & fp2
            c <- length(tmp@bits)

            tmp <- (fp1 | fp2) & !fp2
            a <- length(tmp@bits)

            tmp <- (fp1 | fp2) & !fp1
            b <- length(tmp@bits)

            tmp <- !(fp1 | fp2)
            d <- length(tmp@bits)

            dist <- NULL

            if (method == 'tanimoto') {
              dist <- c / (a+b+c)
            } else if (method == 'euclidean') {
              dist <- sqrt((d+c) / (a+b+c+d))
            } else if (method == 'dice') {
              dist <- c / (.5*a + .5*b + c)
            } else if (method == 'mt') {
              t1 <- c/(size-d)
              t0 <- d/(size-c)
              phat <- ((size-d) + c)/(2*size)
              dist <- (2-phat)*t1/3 + (1+phat)*t0/3
            }

            dist
          })

setGeneric("random.fingerprint",
           function(nbit, on) standardGeneric("random.fingerprint"))
setMethod("random.fingerprint", c("numeric", "numeric"),
          function(nbit, on) {
            if (nbit <= 0) stop("Bit length must be positive integer")
            if (on <= 0) stop("Number of bits to be set to 1 must be positive integer")            
            bits <- sample(1:nbit, size=on)
            new("fingerprint", nbit=nbit, bits=bits, provider="R", folded=FALSE)
          })
