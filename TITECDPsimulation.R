library(Iso)

onetitecdp<-function (PI, target, gam, n, obswin, rate, 
	accrual, minwin, method, cl, scheme, prob) 
{

    	ndose<-length(PI)
    	yvec <- nj <- mj <- ej <- dose.select <- ptox.hat <- rep(0,ndose)
    	len <- length(yvec)
    	userind <- seq(1, len, 1)
    	stop<-0

	x<-2*target
	mu<-target
	u<-0.95
	f<-function(b){
		pbeta(x,mu*b/(1-mu),b)-u
	}
	b0<-uniroot(f,c(0.0001,100))$root
	a0<-mu*b0/(1-mu)

	post.tox<-function(p,a,b,y,m,n,w,delta){
        term=1
	  prior=(((p^(a-1))*(1-p)^(b-1))/(beta(a,b)))
	  for(i in 1:n){
		term=term*(1-w[i]*p)**(1-delta[i])
	  }
	  return(term*(p^y)*((1-p)^(m))*prior);
      }
      
      #the posterior mean of ptox
      posttoxf <- function(p,a,b,y,m,n,w,delta) {p*post.tox(p,a,b,y,m,n,w,delta); }

    	if (accrual == "fixed") {
        next.arrival <- obswin/rate
    	} else if (accrual == "poisson") {
        next.arrival <- rexp(1, rate/obswin)
    	}
        
	  u <- y <- level <- arrival <- weights <- NULL
        cur <- 1
        
        for (i in 1:(n - 1)) {
            arrival <- c(arrival, next.arrival)
            level <- c(level, cur)
         
         	if (method == "uniform") {
            	ynew <- rbinom(1, 1, PI[cur])
            	if (ynew) 
                		unew <- runif(1, 0, obswin)
            	else unew <- Inf
            	y <- c(y, ynew)
            	u <- c(u, unew)
            	utox <- u + arrival
        	}
          	if (method == "exponential") {
                ynew <- rbinom(1, 1, PI[cur])
                if (ynew) {
                  lam <- log(1-PI[cur])/(-obswin)
                        unew <- rexp(1,lam)
				unew <- ifelse(unew>obswin,obswin,unew)
                        }
                else unew <- Inf
                y <- c(y, ynew)
                u <- c(u, unew)
                utox <- u + arrival
            }
        	if (method == "weibull") {
                ynew <- rbinom(1, 1, PI[cur])
                if (ynew) {
                        #unew <- rweibull(1,shape=4,scale=(-(obswin)**4 /log(1-PI[cur]))**0.25)
				alpha.Weibull <- log(log(1 - PI)/log(1 - prob*PI))/log(2)
     				beta.Weibull  <- obswin/(-log(1 - PI))**(1/alpha.Weibull)
      			unew <- min(obswin,rweibull(1,alpha.Weibull[cur],beta.Weibull[cur]))
				unew <- ifelse(unew>obswin,obswin,unew)
                        }
                else unew <- Inf
                y <- c(y, ynew)
                u <- c(u, unew)
                utox <- u + arrival
            }
		if (method == "lognormal")   { 
			ynew <- rbinom(1, 1, PI[cur])
                	if (ynew) {
				sigma <- log(2)/(qnorm(PI) - qnorm(prob * PI))
      			mu    <- (log(1.5) * qnorm(PI) - log(obswin) * qnorm(prob * PI))/(qnorm(PI) - qnorm(prob * PI))
     		 		unew <- rlnorm(1,sigma[cur],mu[cur])
				unew <- ifelse(unew>obswin,obswin,unew)
				}
			else unew <- Inf
                	y <- c(y, ynew)
                	u <- c(u, unew)
                	utox <- u + arrival
            }
		if (method == "gamma") {
			ynew <- rbinom(1, 1, PI[cur])
                	if (ynew) {
	    			f <- function(aa, p = prob, q = 1 - prob) { 
      				qgamma(p, aa, 1)/qgamma(p - p * q, aa, 1) - 2
    				}
    				alpha.Gamma <- beta.Gamma <- rep(0, ndose)
    				for (k in 1:ndose) {
      				alpha.Gamma[k] <- uniroot(f, c(0.01, 10), p = PI[k], q = 1 - prob)$root
      				beta.Gamma[k]  <- qgamma(PI[k], alpha.Gamma[k], 1)/obswin
    					}
				unew <- rgamma(1,alpha.Gamma[cur],beta.Gamma[cur])
				unew <- ifelse(unew>obswin,obswin,unew)
				}
			else unew <- Inf
                	y <- c(y, ynew)
                	u <- c(u, unew)
                	utox <- u + arrival
  		}
        
            if (accrual == "fixed") {
                next.arrival <- next.arrival + obswin/rate
            } else if (accrual == "poisson") {
                next.arrival <- next.arrival + rexp(1, rate/obswin)
            }
            B <- rep(0, length(y))
            B[utox <= next.arrival] <- 1
            censor <- pmin(next.arrival, utox) - arrival
            followup <- pmin(censor, obswin)

		if(scheme=="linear"){
            weights <- followup/obswin
   		}
		if(scheme=="adaptive"){
			support <- sort(followup[y == 1])
			z <- length(support)
			if (z) {
				for (ii in 1:length(followup)) {
      				m <- length(support[support <= followup[ii]])
            			if (!m) 
            				weights[ii] <- followup[ii]/support[1]/(z + 1)
                  		else if (m == z) 
                  			weights[ii] <- (z + (followup[ii] - support[z])/(obswin - 
                      			support[z]))/(z + 1)
                  		else weights[ii] <- (m + (followup[ii] - support[m])/(support[m + 
                  			1] - support[m]))/(z + 1)
               			}
            		} else {
            			weights <- followup/obswin
            		}
			weights[followup==obswin] <- 1
			}

   	 	weights[B == 1] <- 1
    		y1p <- B[1:length(B)]
    		w1p <- weights[1:length(B)]
		
   	 for(i in 1:length(yvec)){
		yvec[i]=sum(y1p[level==i])	
		nj[i]=sum(level==i)
		mj[i] <- sum(w1p[level==i]==1 & y1p[level==i]==0)
		ej[i] <- sum(w1p[level==i]>(minwin/obswin) & y1p[level==i]==0)
	}

	tried=unique(level)
	for(j in 1:length(tried)){
		marginal.tox = integrate(post.tox,lower=0,upper=1,a0,b0,yvec[j],mj[j],nj[j],w1p[level==j],as.numeric(w1p[level==j]==1))$value;
		ptox.hat[j] = integrate(posttoxf,lower=0,upper=1,a0,b0,yvec[j],mj[j],nj[j],w1p[level==j],as.numeric(w1p[level==j]==1))$value/marginal.tox; 
	}

	marginal.tox = integrate(post.tox,lower=0,upper=1,a0,b0,yvec[1],mj[1],nj[1],w1p[level==1],as.numeric(w1p[level==1]==1))$value;
	unsafe1 <-(integrate(post.tox,lower=target,upper=1,a0,b0,yvec[1],mj[1],nj[1],w1p[level==1],as.numeric(w1p[level==1]==1))$value/marginal.tox)>cl
		if(unsafe1){
			stop<-3 
			break
		}
	ptox.hat<-ptox.hat[tried]
	pipost=pava(ptox.hat,w=nj)
		
		lossvec=ifelse(pipost > target, (1-gam)*(pipost-target), gam*(target-pipost))  
		minloss <- min(lossvec)
		havemin <- lossvec==minloss	
		T=lossvec==min(lossvec)
			poss=which(T)
			if(sum(T)==1){
				sugglev=poss
			} else {
				if(all(pipost[poss]>target)){
					sugglev=min(poss)
				} else {
					sugglev=max(poss)
				}
			}
			if(pipost[sugglev]<target & length(tried)<ndose){
				cur<-ifelse(nj[sugglev+1]==0&ej[sugglev]>0,sugglev+1,sugglev)					
			} else {
				cur<-sugglev
			}          
        } #end of for loop


        arrival <- c(arrival, next.arrival)
        level <- c(level, cur)

          	if (method == "uniform") {
            	ynew <- rbinom(1, 1, PI[cur])
            	if (ynew) 
                		unew <- runif(1, 0, obswin)
            	else unew <- Inf
            	y <- c(y, ynew)
            	u <- c(u, unew)
            	utox <- u + arrival
        	}
          	if (method == "exponential") {
                ynew <- rbinom(1, 1, PI[cur])
                if (ynew) {
                  lam <- log(1-PI[cur])/(-obswin)
                        unew <- rexp(1,lam)
				unew <- ifelse(unew>obswin,obswin,unew)
                        }
                else unew <- Inf
                y <- c(y, ynew)
                u <- c(u, unew)
                utox <- u + arrival
            }
        	if (method == "weibull") {
                ynew <- rbinom(1, 1, PI[cur])
                if (ynew) {
                        #unew <- rweibull(1,shape=4,scale=(-(obswin)**4 /log(1-PI[cur]))**0.25)
				alpha.Weibull <- log(log(1 - PI)/log(1 - prob*PI))/log(2)
     				beta.Weibull  <- obswin/(-log(1 - PI))**(1/alpha.Weibull)
      			unew <- min(obswin,rweibull(1,alpha.Weibull[cur],beta.Weibull[cur]))
				unew <- ifelse(unew>obswin,obswin,unew)
                        }
                else unew <- Inf
                y <- c(y, ynew)
                u <- c(u, unew)
                utox <- u + arrival
            }
		if (method == "lognormal")   { 
			ynew <- rbinom(1, 1, PI[cur])
                	if (ynew) {
				sigma <- log(2)/(qnorm(PI) - qnorm(prob * PI))
      			mu    <- (log(1.5) * qnorm(PI) - log(obswin) * qnorm(prob * PI))/(qnorm(PI) - qnorm(prob * PI))
     		 		unew <- rlnorm(1,sigma[cur],mu[cur])
				unew <- ifelse(unew>obswin,obswin,unew)
				}
			else unew <- Inf
                	y <- c(y, ynew)
                	u <- c(u, unew)
                	utox <- u + arrival
            }
       	if (method == "gamma") {
			ynew <- rbinom(1, 1, PI[cur])
                	if (ynew) {
	    			f <- function(aa, p = prob, q = 1 - prob) { 
      				qgamma(p, aa, 1)/qgamma(p - p * q, aa, 1) - 2
    				}
    				alpha.Gamma <- beta.Gamma <- rep(0, ndose)
    				for (k in 1:ndose) {
      				alpha.Gamma[k] <- uniroot(f, c(0.01, 10), p = PI[k], q = 1 - prob)$root
      				beta.Gamma[k]  <- qgamma(PI[k], alpha.Gamma[k], 1)/obswin
    					}
				unew <- rgamma(1,alpha.Gamma[cur],beta.Gamma[cur])
				unew <- ifelse(unew>obswin,obswin,unew)
				}
			else unew <- Inf
                	y <- c(y, ynew)
                	u <- c(u, unew)
                	utox <- u + arrival
  		}

       	weights=rep(1,length(utox))	
    		y1p <- y
    		w1p <- weights
    		 
	for(i in 1:length(yvec)){
		yvec[i]=sum(y1p[level==i])	
		nj[i]=sum(level==i)
		mj[i] <- sum(w1p[level==i]==1 & y1p[level==i]==0)
		ej[i] <- sum(w1p[level==i]>(minwin/obswin) & y1p[level==i]==0)
		}

	tried=unique(level)
	for(j in 1:length(tried)){
		marginal.tox = integrate(post.tox,lower=0,upper=1,a0,b0,yvec[j],mj[j],nj[j],w1p[level==j],as.numeric(w1p[level==j]==1))$value;
		ptox.hat[j] = integrate(posttoxf,lower=0,upper=1,a0,b0,yvec[j],mj[j],nj[j],w1p[level==j],as.numeric(w1p[level==j]==1))$value/marginal.tox; 
	}

	marginal.tox = integrate(post.tox,lower=0,upper=1,a0,b0,yvec[1],mj[1],nj[1],w1p[level==1],as.numeric(w1p[level==1]==1))$value;
	unsafe1 <-(integrate(post.tox,lower=target,upper=1,a0,b0,yvec[1],mj[1],nj[1],w1p[level==1],as.numeric(w1p[level==1]==1))$value/marginal.tox)>cl
		if(unsafe1){
			stop<-3 
		}


	ptox.hat<-ptox.hat[tried]
	pipost=pava(ptox.hat,w=nj)
	
		lossvec=ifelse(pipost > target, (1-gam)*(pipost-target), gam*(target-pipost))  
		minloss <- min(lossvec)
		havemin <- lossvec==minloss	
		T=lossvec==min(lossvec)
			poss=which(T)
			if(sum(T)==1){
				sugglev=poss
			} else {
				if(all(pipost[poss]>target)){
					sugglev=min(poss)
				} else {
					sugglev=max(poss)
				}
			}   
   			mtd=ifelse(stop==3,0,sugglev)
			dose.select[mtd]=dose.select[mtd]+1;    
     return(list(dose.select=dose.select,tox.data=yvec,pt.allocation=nj,duration=max(arrival)+obswin))
}

###Load the function 'titecdp.sim' 
titecdp.sim<-function(ntrial, PI, target, gam, n,obswin, 
			rate, accrual, minwin,method, cl, scheme, prob){
	ndose=length(PI)
	
	d.select<-yo<-no<-matrix(nrow=ntrial,ncol=ndose)
	td<-mtdp<-pcs<-unallowed<-rep(0,ntrial)
	
	for(i in 1:ntrial){
		result<-onetitecdp(PI, target, gam, n, obswin, 
			rate, accrual, minwin,method, cl, scheme, prob) 
		d.select[i,]=result$dose.select
		yo[i,]=result$tox.data
		no[i,]=result$pt.allocation
		pcs[i]=result$dose.select[which.min(abs(PI-target))]
		mtdp[i]=result$pt.allocation[which.min(abs(PI-target))]
		td[i]=result$duration
	}
	cat("Simulation results for Time-to-event CDP design (Wages et al, 2023)\n");
 	cat("targeting a DLT rate of", target,"\n\n");
	cat("True DLT probability:\n")
	cat(round(PI,3),  sep="\t",  "\n")
	cat("MTD Selection percentage:\n")
	cat(formatC(colMeans(d.select)*100, digits=1, format="f"), sep="\t",  "\n")
	cat("Average number of toxicities:\n")
	cat(formatC(colMeans(yo), digits=1, format="f"), sep="\t",   "\n")
	cat("Average number of patients:\n")
	cat(formatC(colMeans(no), digits=1, format="f"), sep="\t",   "\n")
	cat("Percentage of correct selection (PCS):\n")
	cat(round(mean(pcs)*100,2), sep="\t",   "\n")
	cat("Average number of patients treated at mtd:\n")
	cat(round(mean(mtdp),2), sep="\t",   "\n")
	cat("Average trial duration:\n")
	cat(round(mean(td),2), sep="\t",   "\n")
	cat("Percent stopped for safety:\n");
      cat(formatC(100-sum(colMeans(d.select)*100), digits=1, format="f"), sep="\t",   "\n");
		}
##########'titecdp.sim' end here

s1<-c(0.23936641, 0.39043503, 0.55496676, 0.8281236)
s2<-c(0.18943438, 0.31339812, 0.33699722, 0.3623976)
s3<-c(0.09140892, 0.16299757, 0.33395289, 0.5218586)
s4<-c(0.008419951, 0.186018434, 0.37747786, 0.7868993)
s5<-c(0.008471185, 0.09729451, 0.21480590, 0.3349033)
s6<-c(0.112500488, 0.12639304, 0.19461002, 0.3229042)
s7<-c(0.054462330, 0.07246092, 0.08546775, 0.2132631)
s8<-c(0.023416023, 0.066311309, 0.14593070, 0.1838425)
s9<-c(0.58438823, 0.61108281, 0.77066953, 0.7750315)
s10<-c(0.53089357, 0.58016822, 0.68729000, 0.7168440)

target <- 0.2
n <- 24
obswin <- 3
accrual <- "fixed"
rate <- 6
minwin <- 0
ntrial <- 1000
gam <- 0.5
prob <- 0.5
cl <- 0.9
method <- "uniform"
scheme <- "linear"

PI<-s5
set.seed(528957)
titecdp.sim(ntrial, PI, target, gam, n, obswin, rate, accrual, minwin, method, cl, scheme, prob)



