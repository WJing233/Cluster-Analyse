fclustering <- function(x,c=2,noiseDistance=0,m=2,t=1,minchange=0.001,maxsteps=500){
  ### special cases m=1: crisp, 0<m<1: polynomial fuzzifier
  cwn <- c
  if(noiseDistance>0){
    cwn <- c+1
  }  
  w <- matrix(1,nrow=c,ncol=ncol(x))
  if (t>1){
    w <- matrix(1/ncol(x),nrow=c,ncol=ncol(x))
  }
  
  u.old <- matrix(0,nrow=cwn,ncol=nrow(x))
  
  prototypes <- initialise.prototypes(x,c)
  is.converged <- F
  step.count <- 0
  while(!is.converged & step.count<maxsteps){
    step.count <- step.count + 1
    d <- update.distances(x,prototypes,noiseDistance,w,t)
    u <- update.u(d,m)
    is.converged <- max(abs(u-u.old))<minchange
    u.old <- u
    prototypes <- update.prototypes(x,u,noiseDistance,m)
    if (t>1){
      w <- update.weights(x,prototypes,u,m,t)
    } 
  }  
  res <- list(x=x,v=prototypes,u=u,w=w,m=m,t=t)
  #return(step.count)
  return(res)
}

initialise.prototypes <- function(x,c){
  prototypes <- matrix(0,nrow=c,ncol=ncol(x))
  proindex <- sample(1:nrow(x),c,replace=F)
  for (i in 1:nrow(prototypes)){
    for (j in 1:ncol(prototypes)){
      prototypes[i,j] <- x[proindex[i],j]
      if (is.na(prototypes[i,j])){prototypes[i,j] <- mean(x[,j],na.rm=T)}
    }
  } 
  return(prototypes)
}

update.distances <- function(x,prototypes,noiseDistance,w,t){
  nonocl <- 0
  if (noiseDistance>0){
    nonocl <- 1
  }
  d <- matrix(noiseDistance,nrow=nrow(prototypes)+nonocl,ncol=nrow(x))
  for (i in 1:nrow(prototypes)){
    for (j in 1:nrow(x)){
      d[i,j] <- sum(w[i,]^t*(x[j,] - prototypes[i,])^2)
    }
  }
  return(d)
}

g <- function(u,m){
  if (m==1){return(u)}
  if (m>1){return(u^m)}
  return(m*u^2 + (1-m)*u)
}

update.prototypes <- function(x,u,noiseDistance,m){
  nonocl <- 0
  if (noiseDistance>0){
    nonocl <- 1
  }
  prototypes <- matrix(0,nrow=nrow(u)-nonocl,ncol=ncol(x))
  for (i in 1:nrow(prototypes)){
    prototypes[i,] <- t(rep(0,ncol(x)))
    sumu <- 0
    for (j in 1:nrow(x)){
      cij <- g(u[i,j],m)
      sumu <- sumu + cij
      for (k in 1:ncol(prototypes)){
        prototypes[i,k] <- prototypes[i,k] + cij*x[j,k]
      }
    }
    if (sumu>0){
      for (k in 1:ncol(prototypes)){
        prototypes[i,k] <- prototypes[i,k]/sumu
      }
    }
  } 
  return(prototypes)
}

update.u <- function(d,m){
  u <- matrix(0,nrow=nrow(d),ncol=ncol(d))
  if (m==1){
    for (j in 1:(ncol(d))){
      dminj <- d[1,j]
      indmin <- 1
      for (i in 2:(nrow(d))){
        if (dminj>d[i,j]){
          dminj <- d[i,j]
          indmin <- i
        }
      } 
      u[indmin,j] <- 1 
    }    
  } 
  if (m>1){
    for (j in 1:(ncol(d))){
      for (i in 1:(nrow(d))){
        dj <- d[,j]
        if (min(dj)==0){
          dj <- dj + 0.1^30
        }
        u[i,j] <- 1/(sum((d[i,j]/dj)^(1/(m-1))))
        if (is.na(u[i,j])){
          u[i,j] <- 0
        } 
        if (u[i,j]==Inf){
          u[i,j] <- 1
        } 
      }  
    }  
  }
  if (m<1){
    beta <- (1-m)/(1+m) 
    for (j in 1:(ncol(d))){
      dsort <- sort.int(d[,j],index.return=T)$ix
      if (d[dsort[1],j]==0){
        u[dsort[1],j] <- 1
      }else{   
        ul <- (1/beta)-1
        chead <- nrow(d)
        sumdmc <- 0 
        is.below.ul <- F
        while(!is.below.ul & chead>1){
          cv <- (1+(chead-1)*beta)/(d[dsort[chead],j]*sum(1/d[dsort[1:chead],j]) )
          if (cv>beta){
            is.below.ul <- T
          }
          else{
            chead <- chead - 1
          }
        }
        for (k in 1:chead){
          u[dsort[k],j] <- (1/(1-beta))*(((1+(chead-1)*beta)/(d[dsort[k],j]*sum(1/d[dsort[1:chead],j]))) - beta)
          
        }  
      }
    }   
  }
  return(u)
}

update.weights <- function(x,prototypes,u,m,t){
  w <- matrix(0,nrow=nrow(prototypes),ncol=ncol(prototypes))
  expon <- 1/(t-1)
  for (i in 1:nrow(w)){
    for (s in 1:ncol(w)){
      numindenom <- 0
      for (k in 1:(nrow(x))){
        numindenom <- numindenom + g(u[i,k],m)*(x[k,s]-prototypes[i,s])^2
      }
      sumdenom <- 0
      for (r in 1:(ncol(w))){
        sr <- 0
        for (k in 1:(nrow(x))){
          sr <- sr + g(u[i,k],m)*(x[k,r]-prototypes[i,r])^2
        }
        sumdenom <- sumdenom + (numindenom/sr)^expon 
      }
      w[i,s] <- 1/sumdenom
    }
  }
  return(w)
}

plot.clustering.result <- function(x,prototypes,u,dims=c(1,2),pchs=1,xlims=NULL,ylims=NULL,xlab=NULL,ylab=NULL){
  cols <- rep(1,nrow(x))
  for (j in 1:(nrow(x))){
    maxuj <- u[1,j]
    for (i in 2:nrow(u)){
      if (u[i,j]>maxuj){
        maxuj <- u[i,j]
        cols[j] <- i
      }
    }
  }
  xlabt <- xlab
  if (is.null(xlab)){
    xlabt <- paste("Attribute",dims[1])
  }
  
  ylabt <- ylab
  if (is.null(ylab)){
    ylabt <- paste("Attribute",dims[1])
  }
  
  plot(x[,dims[1]],x[,dims[2]],xlim=xlims,ylim=ylims,xlab=xlabt,ylab=ylabt,col=cols,pch=pchs)
  points(prototypes[,dims[1]],prototypes[,dims[2]],col=1:nrow(prototypes),pch=16) 
}


#DDAA

ddaa.one.step <- function(x,prototype,noiseDistance=0,m=2,t=1,minchange=0.001,maxsteps=500){
  ### special cases m=1: crisp, 0<m<1: polynomial fuzzifier
  cwn <- 2
  
  w <- matrix(1,nrow=1,ncol=ncol(x))
  if (t>1){
    w <- matrix(1/ncol(x),nrow=c,ncol=ncol(x))
  }
  
  u.old <- matrix(0,nrow=cwn,ncol=nrow(x))
  
  prototypes <- prototype
  is.converged <- F
  step.count <- 0
  while(!is.converged & step.count<maxsteps){
    step.count <- step.count + 1
    d <- ddaa.update.distances(x,prototypes,noiseDistance,w,t)
    u <- update.u(d,m)
    is.converged <- max(abs(u-u.old))<minchange
    u.old <- u
    prototypes <- ddaa.update.prototype(x,u,m)
    if (t>1){
      w <- ddaa.update.weights(x,prototypes,u,m,t)
    } 
  }  
  res <- list(x=x,v=prototypes,u=u,w=w,m=m,t=t)
  #return(step.count)
  return(res)
}



ddaa.update.distances <- function(x,prototypes,noiseDistance,w,t){
  d <- matrix(noiseDistance,nrow=2,ncol=nrow(x))
  for (j in 1:nrow(x)){
    d[1,j] <- sum(w[1,]^t*(x[j,] - prototypes)^2)
  }
  return(d)
}

ddaa.update.prototype <- function(x,u,m){
  prototype <- rep(0,ncol(x))  
  sumu <- 0
  for (j in 1:nrow(x)){
    cij <- g(u[1,j],m)
    sumu <- sumu + cij
    for (k in 1:length(prototype)){
      prototype[k] <- prototype[k] + cij*x[j,k]
    }
  }
  if (sumu>0){
    for (k in 1:length(prototype)){
      prototype[k] <- prototype[k]/sumu
    }
  }
  
  return(prototype)
}

ddaa.update.weights <- function(x,prototypes,u,m,t){
  w <- matrix(0,nrow=1,ncol=length(prototypes))
  expon <- 1/(t-1)
  
  for (s in 1:ncol(w)){
    numindenom <- 0
    for (k in 1:(nrow(x))){
      numindenom <- numindenom + g(u[1,k],m)*(x[k,s]-prototypes[s])^2
    }
    sumdenom <- 0
    for (r in 1:(ncol(w))){
      sr <- 0
      for (k in 1:(nrow(x))){
        sr <- sr + g(u[1,k],m)*(x[k,r]-prototypes[r])^2
      }
      sumdenom <- sumdenom + (numindenom/sr)^expon 
    }
    w[1,s] <- 1/sumdenom
  }
  
  return(w)
}

cog <- function(x){
  res <- rep(0,ncol(x))  
  for (j in 1:length(res)){
    res[j] <- mean(x[,j]) 
  }
  return(res)
}

multmed <- function(x){
  res <- rep(0,ncol(x))  
  for (j in 1:length(res)){
    res[j] <- median(x[,j]) 
  }
  return(res)
}

plot.ddaa.clustering.result <- function(x,prototypes,u,dims=c(1,2),xlims=NULL,ylims=NULL){
  cols <- rep(2,nrow(x))
  for (j in 1:(nrow(x))){   
    if (u[2,j]>u[1,j]){
      cols[j] <- 1
    }
  }
  plot(x[,dims[1]],x[,dims[2]],xlim=xlims,ylim=ylims,xlab=paste("Attribute",dims[1]),ylab=paste("Attribute",dims[2]),col=cols)
  points(prototypes[dims[1]],prototypes[dims[2]],col=4,pch=16) 
}

plot.ddaa.clustering.result.pca <- function(x,u,xlims=NULL,ylims=NULL,pch=1,label.noise.points=F,shift=0.5,cex.text=1){
  cols <- rep(2,nrow(x))
  for (j in 1:(nrow(x))){   
    if (u[2,j]>u[1,j]){
      cols[j] <- 1
    }
  }
  pca.x <- prcomp(x,center=F,scale=F)
  
  plot(pca.x$x[,1],pca.x$x[,2],xlim=xlims,ylim=ylims,xlab="PC1",ylab="PC2",col=cols,pch=pch)
  
  if (label.noise.points){
    nopos <- subset(1:nrow(x),cols==1)
    text(pca.x$x[nopos,1],pca.x$x[nopos,2]+shift,labels=nopos,cex=cex.text)
  }
}



max.dist.cog <- function(x){
  x.cog <- cog(x)
  max.d <- 0
  for (i in 1:nrow(x)){
    dist.xi <- sum((t(x[i,])-x.cog)^2)
    if (dist.xi>max.d){
      max.d <- dist.xi
    } 	
  }
  return(max.d)
}


init.delta <- function(x,const=5){
  return(const*sqrt(max.dist.cog(x)))
}

entropy.u <- function(pp){
  s <- 0
  for (p in pp){
    if (p>0){
      s <- p*log2(p)
    }
  }
  return(-s)
}

ddaa <- function(x,measures=c("sum.u.real","p.coeff","p.entropy","prototype.change","ambiguous","prototype.displacement"),initial.delta=NULL,const=5,steps=1000,stepwidth=NULL,m=2,t=1,minchange=0.001,maxsteps=500,ambiguous.bounds=c(0.3,0.7),round.sum=T){
  delta0 <- initial.delta
  if (is.null(initial.delta)){
    delta0 <- init.delta(x,const=const)
  }
  
  no.steps <- steps
  if (!is.null(stepwidth)){
    no.steps <- floor(delta0/stepwidth)  
  }
  
  meas.values <- matrix(0,nrow=no.steps,ncol=length(measures))
  colnames(meas.values) <- measures
  
  prototypes <- matrix(0,nrow=no.steps,ncol=ncol(x))
  colnames(prototypes) <- colnames(x)
  
  delta <- rep(0,no.steps)
  
  is.measure <- rep(F,5)
  if (length(subset(measures,measures=="sum.u.real"))>0){
    is.measure[1] <- T
  }
  if (length(subset(measures,measures=="p.coeff"))>0){
    is.measure[2] <- T
  }
  if (length(subset(measures,measures=="p.entropy"))>0){
    is.measure[3] <- T
  }
  if (length(subset(measures,measures=="prototype.change"))>0){
    is.measure[4] <- T
  }
  if (length(subset(measures,measures=="ambiguous"))>0){
    is.measure[5] <- T
  }
  if (length(subset(measures,measures=="prototype.displacement"))>0){
    is.measure[6] <- T
  }
  
  
  p.old <- cog(x)
  p.initial <- p.old
  for (i in 1:no.steps){
    delta[i] <- delta0*(no.steps-i+1)/no.steps
    ddaa1 <- ddaa.one.step(x,p.old,noiseDistance=delta[i]^2,m=m,t=t,minchange=minchange,maxsteps=maxsteps)
    prototypes[i,] <- t(ddaa1$v)
    
    ind.m <- 1
    if (is.measure[1]){
      if (round.sum){
        meas.values[i,ind.m] <- sum(round(t(ddaa1$u[1,])))/nrow(x)
      }else{
        meas.values[i,ind.m] <- sum(t(ddaa1$u[1,]))/nrow(x)
      }
      ind.m <- ind.m+1
    }
    if (is.measure[2]){
      meas.values[i,ind.m] <- sum(ddaa1$u^2)/nrow(x)
      ind.m <- ind.m+1
    }
    if (is.measure[3]){
      meas.values[i,ind.m] <- sum(apply(ddaa1$u,1,entropy.u))/nrow(x)
      ind.m <- ind.m+1
    }
    if (is.measure[4]){
      meas.values[i,ind.m] <- sqrt(sum((p.old-ddaa1$v)^2))
      ind.m <- ind.m+1
    }
    if (is.measure[5]){
      meas.values[i,ind.m] <- sum(ddaa1$u>ambiguous.bounds[1] & ddaa1$u<ambiguous.bounds[2])/(2*nrow(x))
      ind.m <- ind.m+1
    }
    if (is.measure[6]){
      meas.values[i,ind.m] <- sqrt(sum((p.initial-ddaa1$v)^2))
      ind.m <- ind.m+1
    }
    
    p.old <- ddaa1$v
  }
  
  res <- list(v=prototypes,measures=meas.values,delta=delta)
  return(res)
  
}



ddaa.reverse <- function(x,cluster.centre,max.delta,measures=c("sum.u.real","p.coeff","p.entropy","prototype.change","ambiguous","prototype.displacement"),steps=1000,stepwidth=NULL,m=2,t=1,minchange=0.001,maxsteps=500,ambiguous.bounds=c(0.3,0.7)){
  
  no.steps <- steps
  if (!is.null(stepwidth)){
    no.steps <- floor(max.delta/stepwidth)  
  }
  
  meas.values <- matrix(0,nrow=no.steps,ncol=length(measures))
  colnames(meas.values) <- measures
  
  prototypes <- matrix(0,nrow=no.steps,ncol=ncol(x))
  colnames(prototypes) <- colnames(x)
  
  delta <- rep(0,no.steps)
  
  is.measure <- rep(F,5)
  if (length(subset(measures,measures=="sum.u.real"))>0){
    is.measure[1] <- T
  }
  if (length(subset(measures,measures=="p.coeff"))>0){
    is.measure[2] <- T
  }
  if (length(subset(measures,measures=="p.entropy"))>0){
    is.measure[3] <- T
  }
  if (length(subset(measures,measures=="prototype.change"))>0){
    is.measure[4] <- T
  }
  if (length(subset(measures,measures=="ambiguous"))>0){
    is.measure[5] <- T
  }
  if (length(subset(measures,measures=="prototype.displacement"))>0){
    is.measure[6] <- T
  }
  
  
  p.old <- cluster.centre
  p.initial <- p.old
  for (i in 1:no.steps){
    delta[i] <- max.delta*i/no.steps
    ddaa1 <- ddaa.one.step(x,p.old,noiseDistance=delta[i]^2,m=m,t=t,minchange=minchange,maxsteps=maxsteps)
    prototypes[i,] <- t(ddaa1$v)
    
    ind.m <- 1
    if (is.measure[1]){
      meas.values[i,ind.m] <- sum(round(t(ddaa1$u[1,])))/nrow(x)
      ind.m <- ind.m+1
    }
    if (is.measure[2]){
      meas.values[i,ind.m] <- sum(ddaa1$u^2)/nrow(x)
      ind.m <- ind.m+1
    }
    if (is.measure[3]){
      meas.values[i,ind.m] <- sum(apply(ddaa1$u,1,entropy.u))/nrow(x)
      ind.m <- ind.m+1
    }
    if (is.measure[4]){
      meas.values[i,ind.m] <- sqrt(sum((p.old-ddaa1$v)^2))
      ind.m <- ind.m+1
    }
    if (is.measure[5]){
      meas.values[i,ind.m] <- sum(ddaa1$u>ambiguous.bounds[1] & ddaa1$u<ambiguous.bounds[2])/(2*nrow(x))
      ind.m <- ind.m+1
    }
    if (is.measure[6]){
      meas.values[i,ind.m] <- sqrt(sum((p.initial-ddaa1$v)^2))
      ind.m <- ind.m+1
    }
    
    p.old <- ddaa1$v
  }
  
  res <- list(v=prototypes,measures=meas.values,delta=delta)
  return(res)
  
}



plot.measures <- function(vmd,measure.names=c("Sum","PC","PE","Jumps","Amb","Displ"),dSum=T,cols=1:7,lwd=2,lty=1:7,legend.loc=NULL){
  plot(NULL,xlim=c(min(vmd$delta),max(vmd$delta)),ylim=c(0,1),xlab=expression("delta"),ylab="",main="")
  nmwds <- length(measure.names)
  for (i in 1:nmwds){
    const <- 1
    if (measure.names[i]=="Jumps" | measure.names[i]=="PE" |  measure.names[i]=="Displ"){
      const <- max(vmd$measures[,i])
    }
    points(vmd$delta,vmd$measures[,i]/const,type="l",col=cols[i],lwd=lwd,lty=lty[i])
  }
  if (dSum){
    ds <-  c(1,vmd$measures[1:(nrow(vmd$measures)-1)]) - vmd$measures[,1] 
    ds <- ds/max(ds)
    points(vmd$delta,ds,type="l",col=cols[length(cols)],lwd=lwd,lty=lty[length(lty)])
  }
  mn <- measure.names
  if (dSum){
    mn <- c(mn,"dSum")
  }
  if (!is.null(legend.loc)){
    if (length(legend.loc)==1){
      legend(legend.loc,legend=mn,col=cols,lwd=lwd,lty=lty)
    }else{
      legend(legend.loc[1],legend.loc[2],legend=mn,col=cols,lwd=lwd,lty=lty)
    }
  }
}


plot.ddaa.result <- function(x,prototypes,dims=c(1,2),xlims=NULL,ylims=NULL,cols=1,prototype.colour=2,xlab=NULL,ylab=NULL,add=F){
  if (!add){
    xlabt <- xlab
    if (is.null(xlab)){
      xlabt <- paste("Attribute",dims[1])
    }
    
    ylabt <- ylab
    if (is.null(ylab)){
      ylabt <- paste("Attribute",dims[1])
    }
    plot(x[,dims[1]],x[,dims[2]],xlim=xlims,ylim=ylims,xlab=xlabt,ylab=ylabt,col=cols)
  }	
  points(prototypes[1,dims[1]],prototypes[1,dims[2]],col=prototype.colour,pch=15) 
  points(prototypes[nrow(prototypes),dims[1]],prototypes[nrow(prototypes),dims[2]],col=prototype.colour,pch=17) 
  points(prototypes[,dims[1]],prototypes[,dims[2]],col=prototype.colour,type="l")   
}

sele<-function(x,kk,colnumber,rownumber){
  if(x==1){
    ss<-cog(iris.num[1:rownumber,1:colnumber])
    return(ss)
  }else if(x==2){
    ss<-iris.num[kk,]
    return(ss)
  }
}








library(shiny)
library(shinydashboard)
library(ggplot2)


ui <- dashboardPage(
  dashboardHeader(title="Cluster Analysis"),
  dashboardSidebar(hr(), 
                   sidebarMenu(
                     menuItem(text = "Settings",tabName = "settings",icon = icon("cog")),
                     menuItem(text = "Import",tabName = "import",icon = icon("table")),
                     menuItem(text ="Plot",tabName = "plot",icon = icon("chart-line")),
                     menuItem(text = "Download",tabName = "download",icon=icon("download")),
                     menuItem(text = "About",tabName = "about",icon=icon("clipboard")),
                     disable = FALSE   
                   ),
                   hr()
  ),
  
  dashboardBody(
    tabItems(
      tabItem("settings",
              fluidRow(
                box(width = 12,
                    solidHeader = TRUE,
                    title = "Settings",
                    h3("Dataset"),
                    textInput("numberofcolumn",
                              "number of column",
                              value = "4"
                    ),
                    
                    hr(),
                    
                    fluidRow(
                      column(width = 12,
                             box(width = NULL, 
                                 status = "warning",
                                 collapsible = TRUE,
                                 title = "Settings of Plot1",
                                 solidHeader = FALSE,
                                 textInput("xlabel",
                                           "enter the x label",
                                 ),
                                 
                                 textInput("ylabel",
                                           "enter the y label",
                                 ),
                                 
                             )
                      )
                    ),
                    
                    fluidRow(
                      column(
                        width=12,
                        box(
                          width=NULL,
                          status = "warning",
                          collapsible = TRUE,
                          title="Settings of Plot4",
                          solidHeader = FALSE,
                          textInput("start",
                                    "enter the start Row of the cluster that you selected",
                                    value = "1"
                          ),
                          
                          textInput("end",
                                    "enter the end Row of the cluster that you selected",
                                    value = "50"
                          ),
                        )
                      ),
                    ),
                    
                    fluidRow(
                      column(
                        width=12,
                        box(
                          width=NULL,
                          status = "warning",
                          collapsible = TRUE,
                          title="Settings of Plot6",
                          solidHeader = FALSE,
                          radioButtons("clucen","Choose the type of the startpoint",
                                       choices =list("average"=1,"new start point"=2),
                                       selected = 1 ),
                          textInput("clusterrownumber","please enter the number of the sample from dataset",
                                    placeholder = ""),
                          
                        )
                      ),
                      
                    ),
                    
                    
                )
              )
      ),
      
      
      tabItem("import",
              fluidRow(
                box(width = 12,
                    solidHeader = TRUE,
                    title = "Import",
                    fileInput("filename","File upload:",accept=c(".csv")),
                    br(),
                    
                ),
                box(width=12,
                    collapsible = TRUE,
                    title="Details",
                    solidHeader = TRUE,
                    tableOutput("table1")
                ),
              ),
              
              
              
              
      ),
      
      
      
      tabItem("plot",
              fluidRow(
                
                tabBox(width = 12,
                       tabPanel("Plot1",
                                "This shows the Plot of FCM result",
                                hr(),
                                tabPanel("",plotOutput("image1")),
                                solidHeader=T,
                                width = 12
                       ),
                       
                       
                       
                       
                       tabPanel("Plot2",
                                "This shows the Plot of DDAA result for the dataset",
                                hr(),
                                tabPanel("",plotOutput("image2")),
                                solidHeader=T,
                                width = 12
                       ),
                       
                       tabPanel("Plot3",
                                "DDAA Prototype movement",
                                hr(),
                                tabPanel("",plotOutput("image3")),
                                solidHeader=T,
                                width = 12
                       ),
                       
                       tabPanel("Plot4",
                                "Validating single clusters",
                                hr(),
                                tabPanel("",plotOutput("image5")),
                                solidHeader=T,
                                width = 12 
                       ),
                       
                       tabPanel("Plot5",
                                "starting the cluster from the middle",
                                hr(),
                                tabPanel("",plotOutput("image7")),
                                solidHeader=T,
                                width = 12 
                       ),
                       
                       tabPanel("Plot6",
                                "prototypemovements",
                                hr(),
                                tabPanel("",plotOutput("image6")),
                                solidHeader=T,
                                width = 12 
                       )
                ),
                
                
              ),
              
              
              
      ),
      
      tabItem("about",
              fluidRow(
                
                box(
                  width = 12,
                  solidHeader = TRUE,
                  title = "About",
                  h3("This is a application that can show the situation about the dataset and aimed for the cluster analysis."),
                  hr(),
                  box(
                    width = 12,
                    title = "Step 1",
                    "Set up the name of the x-y axis and enter the target clustering. It can also set up the startpoint of the clustering",
                    
                  ),
                  hr(),
                  box(
                    width = 12,
                    title = "Step 2",
                    "Upload the csv data of the clustering and the details will be showed under it",
                    
                  ),
                  hr(),
                  box(
                    width = 12,
                    title = "Step 3",
                    "The result will in the form of Plots displayed through the caculating",
                    
                  ),
                  
                )
              )
      ),
      
      tabItem("download",
              fluidRow(
                box(
                  width = 12,
                  solidHeader = TRUE,
                  title = "Download",
                  column(width = 6,
                         downloadButton(outputId = "downpic1", label = "Download the plot1"),
                  ),
                  column(width = 6,
                         downloadButton(outputId = "downpic2", label = "Download the plot2"),
                  ),
                  
                  column(width = 12,
                         hr(),
                  ),
                  
                  column(width = 6,
                         downloadButton(outputId = "downpic3", label = "Download the plot3"),
                  ),
                  column(width = 6,
                         downloadButton(outputId = "downpic4", label = "Download the plot4"),
                  ),
                  
                  column(width = 12,
                         hr(),
                  ),
                  
                  column(width = 6,
                         downloadButton(outputId = "downpic5", label = "Download the plot5"),
                  ),
                  column(width = 6,
                         downloadButton(outputId = "downpic6", label = "Download the plot6"),
                  ),
                  
                  column(width = 12,
                         hr(),
                  ),
                  
                )
              )
      )
      
      
      
    )
    
    
  )
  
  
  
)

server <- function(input,output){
  
  drawpic1<-reactive({
    if(!is.null(input$filename)){
      file_spec <- input$filename
      bb <- read.csv(file_spec$datapath,header = TRUE,sep=",")
      iris.num <<- bb[,1:input$numberofcolumn]
      
      fcm.iris <<- fclustering(iris.num,c=3,noiseDistance=0,m=2,t=1,minchange=0.001,maxsteps=500)
      
      pchs <<- rep(21,nrow(iris.num))
      for (j in 1:(nrow(iris.num))){
        maxuj <- fcm.iris$u[1,j]
        for (i in 2:nrow(fcm.iris$u)){
          if (fcm.iris$u[i,j]>maxuj){
            maxuj <- fcm.iris$u[i,j]
            if (i==2){
              pchs[j] <- 22
            }else{
              pchs[j] <- 24
            }
          }
        }
      }
      
    }
    
  })
  
  output$image1<-renderPlot({
    if(is.null(data_file())){
      return()
    }
    else{
      drawpic1()
      plot.clustering.result(iris.num,fcm.iris$v,fcm.iris$u,dims=c(3,4),pchs=pchs,xlab=input$xlabel,ylab=input$ylabel)
      
    }
    
  })
  
  
  drawpic2<-reactive({
    if(is.null(input$filename)){
      return()
    }
    else{
      file_spec <- input$filename
      bb <- read.csv(file_spec$datapath,header = TRUE,sep=",")
      diris <<- ddaa(iris.num,const=1)
      
    }
  })
  
  output$image2<-renderPlot({
    if(is.null(data_file())){
      return()
    }
    else{
      drawpic2()
      plot.measures(diris,dSum=F,legend.loc=c(2.9,0.33),lwd=3)
      
    }
    
  })
  
  
  drawpic3<-reactive({
    if(is.null(input$filename)){
      return()
    }
    else{
      file_spec <- input$filename
      bb <- read.csv(file_spec$datapath,header = TRUE,sep=",")
      
    }
  })
  
  output$image3<-renderPlot({
    if(is.null(data_file())){
      return()
    }
    else{
      drawpic3()
      plot.ddaa.result(iris.num[,1:input$numberofcolumn],rbind(cog(iris.num[,1:input$numberofcolumn]),diris$v),dims=c(3,4),xlab=input$xlabel,ylab=input$ylabel)
      
    }
    
    
  })
  
  
  drawpic5<-reactive({
    if(is.null(input$filename)){
      return()
    }
    else{
      file_spec <- input$filename
      bb <- read.csv(file_spec$datapath,header = TRUE,sep=",")
      
      setosa.max.dist <- max.dist.cog(iris.num[input$start:input$end,])

      setosa.cog <- cog(iris.num[input$start:input$end,])
      
      iris.seto <<- ddaa.reverse(iris.num,setosa.cog,2*setosa.max.dist,steps=100)
      
    }
  })
  
  output$image5<-renderPlot({
    if(is.null(data_file())){
      return()
    }
    else{
      drawpic5()
      plot.measures(iris.seto,dSum=F,legend.loc=c(2.55,0.35),lwd=3)
      
    }
    
  })
  
  drawpic7<-reactive({
    if(is.null(input$filename)){
      return()
    }
    else{
      file_spec <- input$filename
      bb <- read.csv(file_spec$datapath,header = TRUE,sep=",")
      
      c2.max.dist <- max.dist.cog(iris.num[is.clu2,])
      #2.670372
      
      c2.cog <- cog(iris.num[is.clu2,])
      
      iris.c2 <<- ddaa.reverse(iris.num,c2.cog,2*c2.max.dist,steps=100)
      
    }
  })
  
  output$image7<-renderPlot({
    if(is.null(data_file())){
      return()
    }
    else{
      drawpic7()
      plot.measures(iris.c2,dSum=F,legend.loc=c(4.4,0.7),lwd=3)      
      
    }
    
  })
  
  
  drawpic6<-reactive({
    if(is.null(input$filename)){
      return()
    }
    else{
      cluscent<<-sele(input$clucen,input$clusterrownumber,input$numberofcolumn,input$end)
      
      file_spec <- input$filename
      bb <- read.csv(file_spec$datapath,header = TRUE,sep=",")
      
      
      
    }
  })
  
  output$image6<-renderPlot({
    if(is.null(data_file())){
      return()
    }
    else{
      drawpic6()
      plot.ddaa.result(iris.num[,1:input$numberofcolumn],rbind(cluscent,iris.seto$v),dims=c(3,4),xlab=input$xlabel,ylab=input$ylabel)
      plot.ddaa.result(iris.num[,1:input$numberofcolumn],rbind(cog(iris.num[is.clu2,]),iris.c2$v),dims=c(3,4),prototype.colour=4,add=T,)
      
    }
    
  })
  
  
  data_file<-reactive({
    if(is.null(input$filename)){
      return()
    }
    else{
      file_spec<-input$filename
      aa<-read.table(file_spec$datapath,header = TRUE,sep=",")
      return(aa)
    }
  })
  
  output$table1<-renderTable({
    if(is.null(data_file())){
      return()
    }
    else{
      is.clu2 <- rep(F,ncol(fcm.iris$u))
      for (j in 1:(ncol(fcm.iris$u))){
        clu.ind.j <- 1
        maxuj <- fcm.iris$u[1,j]
        for (i in 2:nrow(fcm.iris$u)){
          if (fcm.iris$u[i,j]>maxuj){
            maxuj <- fcm.iris$u[i,j]
            clu.ind.j <- i
          }
        }
        is.clu2[j] <- clu.ind.j==2
      }
      
      data_file()
    }  
    
  })
  
  output$downpic1<- downloadHandler(filename = function(){
    paste("plot1","png",sep = ".")
  },
  content = function(file){
    png(file)
    plot.clustering.result(iris.num,fcm.iris$v,fcm.iris$u,dims=c(3,4),pchs=pchs,xlab=input$xlabel,ylab=input$ylabel)
    dev.off()
  },
  contentType = 'png'
  )
  
  
  output$downpic2<- downloadHandler(filename = function(){
    paste("plot2","png",sep = ".")
  },
  content = function(file){
    png(file)
    plot.measures(diris,dSum=F,legend.loc=c(2.9,0.33),lwd=3)
    dev.off()
  },
  contentType = 'png'
  )
  
  
  output$downpic3<- downloadHandler(filename = function(){
    paste("plot3","png",sep = ".")
  },
  content = function(file){
    png(file)
    plot.ddaa.result(iris.num[,1:input$numberofcolumn],rbind(cog(iris.num[,1:input$numberofcolumn]),diris$v),dims=c(3,4),xlab=input$xlabel,ylab=input$ylabel)
    dev.off()
  },
  contentType = 'png'
  )
  
  output$downpic4<- downloadHandler(filename = function(){
    paste("plot4","png",sep = ".")
  },
  content = function(file){
    png(file)
    plot.measures(iris.seto,dSum=F,legend.loc=c(2.55,0.35),lwd=3)
    dev.off()
  },
  contentType = 'png'
  )
  
  output$downpic5<- downloadHandler(filename = function(){
    paste("plot5","png",sep = ".")
  },
  content = function(file){
    png(file)
    plot.measures(iris.c2,dSum=F,legend.loc=c(4.4,0.7),lwd=3)
    dev.off()
  },
  contentType = 'png'
  )
  
  output$downpic6<- downloadHandler(filename = function(){
    paste("plot6","png",sep = ".")
  },
  content = function(file){
    png(file)
    plot.ddaa.result(iris.num[,1:input$numberofcolumn],rbind(cog(iris.num[1:input$end,1:input$numberofcolumn]),iris.seto$v),dims=c(3,4),xlab=input$xlabel,ylab=input$ylabel)
    plot.ddaa.result(iris.num[,1:input$numberofcolumn],rbind(cog(iris.num[is.clu2,]),iris.c2$v),dims=c(3,4),prototype.colour=4,add=T,)
    dev.off()
  },
  contentType = 'png'
  )
  
}





# Run the application 
shinyApp(ui = ui, server = server)
