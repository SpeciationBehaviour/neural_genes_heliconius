## Figure 1 & Supplementary Figure 1 - PLOTS ##
Alexander E. Hausmann (alexander_hausmann@gmx.net), March 2020

## Initial Setup

Setting a seed (basically only relevant for reproducing specific jitter pattern in the plots).
```{r}
set.seed(42)
```

Set working directory (change to respective device)
```{r}
setwd("C:/Users/Hausmann/Desktop/Rossi_et_al_2020/")
```

For reading PNGs
```{r}
suppressMessages(suppressWarnings(library(png)))
```

For displaying PNGs
```{r}
suppressMessages(suppressWarnings(library(grid)))
```

Read in data created in analysis script
```{r}
load("analyses_fig1_suppl_fig_1.RData")
```

Read in *cydno* and *melpomene* photo
```{r}
cydno<-readPNG("cydno.png")
melpomene<-readPNG("melpomene.png")
```

## Define function for ternary plots

In the following, the function `ternary_choice` will be defined, which later will be called multiple times to produce the plots.

First, some subfunctions will be defined which later will be called within `ternary_choice`.

The left (skewed) ternary axis will be further referred to as x1, the bottom one as x2 and the right (skewed) one as x3. All functions or variables relating to one of these axes will carry the respective abbreviation in their name.

### Subfunction definitions

Tick drawing functions

```{r}
tickx1<-function(xcoordo,ticks_inner_outer,tick_length,spaco2,...){
  #calculate where diagonals cut throw edge line
  x1<-((sqrt(3)*(1-xcoordo)-2*spaco2)/(2*sqrt(3)))-0.5
  #plug into one of the lines
  y1<-sqrt(3)*(x1+0.5)+2*spaco2
  x2<-((sqrt(3)*(1-xcoordo)-(2*(spaco2+tick_length*sin(60*(pi/180)))))/(2*sqrt(3)))-0.5
  y2<-sqrt(3)*(x2+0.5)+2*(spaco2+tick_length*sin(60*(pi/180)))
  segments(x1,y1,x2,y2,...)
}

tickx2<-function(xcoordo,ticks_inner_outer,tick_length,spaco2,...){
  if(ticks_inner_outer=="outer"){
    add_space<-spaco2
  } else{
    add_space<-0
  }
  xtransl<-xcoordo-add_space*(3^(-(1/2)))-0.5
  ynew<-(-add_space-tick_length*sin(60*(pi/180)))
  xnew<-xcoordo+ynew*(3^(-(1/2)))-0.5
  segments(xtransl,-add_space,xnew,ynew,...)
}

tickx3<-function(xcoordo,ticks_inner_outer,tick_length,spaco2,...){
  if(ticks_inner_outer=="outer"){
    add_space<-spaco2
  } else{
    add_space<-0
  }
  xtransl<-1-(xcoordo/2)+(add_space/sin(60*(pi/180)))-0.5
  yboth<-xcoordo*(sqrt(3)/2)
  xnew<-xtransl+tick_length
  segments(xtransl,yboth,xnew,yboth,...)
}
```

Axis labels functions

```{r}
axx1<-function(xcoordo,ticks_inner_outer,tick_length,axis_d,spaco2,
               y_over_x_ext1,...){
  x2<-((sqrt(3)*(1-xcoordo)-(2*(spaco2+(tick_length+axis_d)*sin(60*(pi/180)))))/(2*sqrt(3)))-0.5
  y2<-sqrt(3)*(x2+0.5)+2*(spaco2+(tick_length+axis_d)*sin(60*(pi/180)))
  text(x2,y2,xcoordo,srt=360-(atan(sqrt(3)*y_over_x_ext1)*180)/pi,...)
}

axx2<-function(xcoordo,ticks_inner_outer,tick_length,axis_d,spaco2,
               y_over_x_ext1,...){
  if(ticks_inner_outer=="outer"){
    add_space<-spaco2
  } else{
    add_space<-0
  }
  ynew<-(-add_space-(tick_length+axis_d)*sin(60*(pi/180)))
  xnew<-xcoordo+ynew*(3^(-(1/2)))-0.5
  text(xnew,ynew,xcoordo,srt=(atan(sqrt(3)*y_over_x_ext1)*180)/pi,...)
}

axx3<-function(xcoordo,ticks_inner_outer,tick_length,axis_d,spaco2,...){
  if(ticks_inner_outer=="outer"){
    add_space<-spaco2
  } else{
    add_space<-0
  }
  xtransl<-1-(xcoordo/2)+(add_space/sin(60*(pi/180)))-0.5
  xnew<-xtransl+(tick_length+axis_d)
  ynew<-xcoordo*(sqrt(3)/2)
  text(xnew,ynew,xcoordo,...)
}
```

Coordinate system grid functions

```{r}
ax1<-function(xcoordo,...){
  xstart<-1-xcoordo-0.5
  ystart<-0
  xnew<-(0.5-xcoordo/2)-0.5
  ynew<-(1-xcoordo)*(sqrt(3)/2)
  segments(xstart,ystart,xnew,ynew,...)
}

ax2<-function(xcoordo,...){
  ynew<-(1-xcoordo)*(sqrt(3)/2)
  xnew<-((ynew+sqrt(3)*xcoordo)/sqrt(3))-0.5
  segments(xcoordo-0.5,0,xnew,ynew,...)
}

ax3<-function(xcoordo,...){
  xstart<-(1-xcoordo/2)-0.5
  ystart<-xcoordo*(sqrt(3)/2)
  xnew<-(-abs(xstart))
  ynew<-xcoordo*(sqrt(3)/2)
  segments(xstart,ystart,xnew,ynew,...)
}
```

Transfer dots to ternary space

```{r}
transf_x<-function(left,bottom,right){
  return(0.5*((2*bottom+right)/(left+bottom+right))-0.5)
}
transf_y<-function(left,bottom,right){
  return((sqrt(3)/2)*(right/((left+bottom+right))))
}
```

### Main function definition

```{r}
ternary_choice<-function(
  ### Output file path. If "", no png will be made. If path is given, 
  ### resolution and extent of png will be optimized automatically
  output_file="",
  ### Pass data frame with three columns: 
  #First: left axis, Second: bottom axis, Third: right axis
  #If c(), all display of raw data skipped.
  data_fr,
  ### Pass vector of colours for dots
  dot_cols,
  ### Vector of dot sizes
  dot_sizes,
  ### Transparency of dots
  transp=1/3,
  ### Jitter applied to dots
  jitter_dots=0,
  ### Estimators and CrIs. This is a list. One entry of the
  #list stands for one species and gives estimators in all directions.
  #One entry of the list is a data frame of three rows and three columns:
  #First column: estimators, second column: lower boundary, third column: upper boundary
  #First row: towards left axis, second row: towards bottom axis, third row: towards right axis
  #If you don't want this feature, just put to c()
  est_coord,
  ### Colours used for estimator and CrI hexagon. A vector of equal length as est_coord
  est_cols="red",
  ### Helping thin lines from estimator to axes
  est_line=T,
  ### Whether to add contours
  cont_yn=T,
  ### List of data points for each subset of contours.
  ### E.g. rows 1,2,3 should get one contour and rows 4,5,6 another. Then set this to
  ### list(c(1,2,3),c(4,5,6))
  cont_list,
  ### Weights for contours. Give as a list of same format as cont_list. Only use integers.
  ### This increases the "number" of data points and therefore the
  ### logic of the density estimation. Density is then estimated from your "raw"
  ### data and not the summarized table (your summarized table gets extended by
  ### weights in this step).
  cont_wts,
  ### nlevels for contours
  nlevels_cont=10,
  ### Vector of colours for each contour set
  cont_cols=c("black"),
  ### Whether to remove all parts of the dots hanging over the edges of the triangle
  remove_outer=F,
  ### Colour of background
  bg_col=adjustcolor("red",0),
  ### Whether ticks go directly on inner triangle or on outer
  ticks_inner_outer="outer",
  ### Space distance of outer triangle
  spaco2=0,
  ### Distance of arrows from triangle
  disto=0.15,
  ### Half the length of an arrow
  half_length_arrow=0.3,
  ### Arrow line width
  arrow_lwd=1,
  ### Distance between axis title and arrow
  disto_lab=0.03,
  ### Length of ticks
  tick_length1=0.03,
  ### Space around triangle
  space_around=c(0.25,0.2,0.1,0.2),
  ### Positions at which helping lines are drawn in triangle
  line_at=seq(0,1,0.1),
  ### Positions at which ticks are drawn
  tck_at=seq(0,1,0.1),
  ### Axis values and their position
  axis_l=seq(0,1,0.1),
  ### Size of axis labels
  axis_cex=0.7,
  ### Distance between axis label and tick end
  axis_d1=0.025,
  ### Names of three axis titles, from left, to bottom, to right
  title_names=c("sp1","sp2","sp3"),
  ### Size of titles
  title_cex=1,
  ### Vector naming scales to plot. 1 stands for scale from left axis to
  ### right bottom corner, 2 for axis from bottom to top, 3 for axis 
  ### from right axis to left bottom corner.
  scales_to_plot=1:3
){
  
  #######################################################
  
  #Open png if wished. Set dimensions such that the angles of the triangle are 60 degrees.
  if(output_file!=""){
    png(output_file,width=1500,
        height=((((sqrt(3))/2)+sum(space_around[c(1,3)]))/(1+sum(space_around[c(2,4)])))*1500,
        res=300)
  }
  
  #Replace if some arguments were not passed by user.
  if(missing(dot_cols)&length(data_fr)>0){
    dot_cols<-rep("gray",nrow(data_fr))
  }
  if(missing(dot_sizes)&length(data_fr)>0){
    dot_sizes<-rep(1,nrow(data_fr))
  }
  if(missing(cont_list)&length(data_fr)>0){
    cont_list<-list(1:nrow(data_fr))
  }
  if(missing(cont_wts)&cont_yn==T){
    cont_wts<-lapply(1:length(cont_list),function(x) rep(1,length(cont_list[[x]])))
  }
  
  #Calculate ternary position of points (and jitter, if wished)
  if(length(data_fr)>0){
    new_x<-jitter(transf_x(left=data_fr[,1],bottom=data_fr[,2],right=data_fr[,3]),jitter_dots)
    new_y<-jitter(transf_y(left=data_fr[,1],bottom=data_fr[,2],right=data_fr[,3]),jitter_dots)
  }
  
  #Remove outer white space
  #Safe first old mar setting
  old_mar<-par("mar")
  par(mar=c(0,0,0,0))
  
  #Open plot
  plot(1,1,xlim=c(-space_around[2]-0.5,1+space_around[4]-0.5),
       ylim=c(-space_around[1],((sqrt(3))/2)+space_around[3]),xaxs="i",yaxs="i",
       xaxt="n",yaxt="n",xlab="",ylab="",type="n",bty="n")
  
  #This helps later with positioning. 
  #It's the ratio of width and height of plotting device times
  #ratio of width and height of Carthesian coordinate system.
  y_over_x_ext<-(par("pin")[2]/par("pin")[1])*(diff(par("usr")[1:2])/diff(par("usr")[3:4]))
  
  #Add desired background colour
  rect(-1000,-1000,1000,1000,border=NA,col=bg_col)
  
  #Add white triangle
  polygon(c(0-0.5,1-0.5,0.5-0.5),c(0,0,sqrt(3)/2),col="white",border=NA)
  
  #The coordinate system grid
  ax1(line_at,col="grey70",lwd=0.5)
  ax2(line_at,col="grey70",lwd=0.5)
  ax3(line_at,col="grey70",lwd=0.5)
  
  #The scales for each axis
  if(1%in%scales_to_plot){
    segments(1-0.5,0,0.25-0.5,0.5*(sqrt(3)/2),lty="dashed")
  }
  if(2%in%scales_to_plot){
    segments(0.5-0.5,0,0.5-0.5,sqrt(3)/2,lty="dashed")
  }
  if(3%in%scales_to_plot){
    segments(0-0.5,0,0.75-0.5,0.5*(sqrt(3)/2),lty="dashed")
  }
  # (REMOVED: These would be the 1/3 lines)
  # ax1(1/3,col="red",lty="dashed")
  # ax2(1/3,col="blue",lty="dashed")
  # ax3(1/3,col="green",lty="dashed")
  
  #If outer area should be removed
  if(remove_outer){
    #Add points
    if(length(data_fr)>0){
      points(new_x,new_y,pch=21,
             bg=adjustcolor(dot_cols,transp),cex=dot_sizes)
    }
    #Remove all outside triangle
    polygon(c(-1000,1000,1000,-1000),
            c(-1000,-1000,-spaco2,-spaco2),col=bg_col,border=NA)
    m_new<-(-sqrt(3))
    t_new<-(-spaco2)-(m_new*(1+spaco2*(3/sqrt(3))))
    polygon(c(1+spaco2*(3/sqrt(3)),1000,1000+(1000-t_new)/m_new,(1000-t_new)/m_new)-0.5,
            c(-spaco2,-spaco2,1000,1000),col=bg_col,border=NA)
    t_new2<-(-spaco2)-(-m_new*(-spaco2*(3/sqrt(3))))
    polygon(c(-1000,-spaco2*(3/sqrt(3)),(1000-t_new2)/(-m_new),-1000+spaco2*(3/sqrt(3))+(1000-t_new2)/(-m_new))-0.5,
            c(-spaco2,-spaco2,1000,1000),col=bg_col,border=NA)
  }
  
  #Draw "axes"
  segments(c(0,0.5,0)-0.5,
           c(0,sqrt(3)/2,0),
           c(0.5,1,1)-0.5,
           c(sqrt(3)/2,0,0),"black",lwd=0.75)
  
  #Draw outer triangle
  #get left bottom corner
  #y=0*x-spaco2
  #y=(sqrt(3)/4)/0.75
  #Intersection: x=-spaco2*(3/sqrt(3))
  #Same logic for bottom right
  segments(c(-spaco2*(3/sqrt(3)),0.5,-spaco2*(3/sqrt(3)))-0.5,
           c(-spaco2,(sqrt(3)/2)+2*spaco2,-spaco2),
           c(0.5,1+spaco2*(3/sqrt(3)),1+spaco2*(3/sqrt(3)))-0.5,
           c((sqrt(3)/2)+2*spaco2,-spaco2,-spaco2),lwd=1.5)
  
  #Add ticks
  tickx1(tck_at,ticks_inner_outer=ticks_inner_outer,col="black",lwd=1,
         tick_length=tick_length1,spaco2=spaco2)
  tickx2(tck_at,ticks_inner_outer=ticks_inner_outer,col="black",lwd=1,
         tick_length=tick_length1,spaco2=spaco2)
  tickx3(tck_at,ticks_inner_outer=ticks_inner_outer,col="black",lwd=1,
         tick_length=tick_length1,spaco2=spaco2)
  
  #Add axis labels
  if(axis_cex>0){
    axx1(xcoordo=axis_l,ticks_inner_outer=ticks_inner_outer,
         tick_length=tick_length1,axis_d=axis_d1,cex=axis_cex,spaco2=spaco2,
         y_over_x_ext1=y_over_x_ext)
    axx2(xcoordo=axis_l,ticks_inner_outer=ticks_inner_outer,
         tick_length=tick_length1,axis_d=axis_d1,cex=axis_cex,spaco2=spaco2,
         y_over_x_ext1=y_over_x_ext)
    axx3(xcoordo=axis_l,ticks_inner_outer=ticks_inner_outer,
         tick_length=tick_length1,axis_d=axis_d1,cex=axis_cex,spaco2=spaco2)
  }
  
  #Add points (if not added before already)
  if(!remove_outer){
    if(length(data_fr)>0){
      points(new_x,new_y,pch=21,
             bg=adjustcolor(dot_cols,transp),cex=dot_sizes)
    }
  }
  
  #Add contours, if wished.
  if(cont_yn){
    if(length(data_fr)>0){
      #Call Ternary package and set ternDirection to 1.
      if(!"Ternary"%in% rownames(installed.packages())){
        install.packages("Ternary")
      }
      suppressMessages(suppressWarnings(library(Ternary)))
      options(ternDirection=1)
      #Go one by one through groups
      for(conturo in 1:length(cont_list)){
        cont_sub<-data_fr[cont_list[[conturo]],c(3,2,1)]
        #Inflate table by weights
        cont_sub_new<-cont_sub[-(1:nrow(cont_sub)),]
        for(multipl in 1:nrow(cont_sub)){
          for(x in 1:cont_wts[[conturo]][multipl]){
            cont_sub_new<-rbind(cont_sub_new,cont_sub[multipl,])
          }
        }
        #Calculate Ternary density (using the weighted/bloated data set)
        TernaryDensityContour(coordinates=cont_sub_new,
                              drawlabels=F,col=cont_cols[conturo],direction=1,
                              nlevels=nlevels_cont)
      }
    }
  }
  
  #Positions of text
  x_textx1<-0.25-sin(60*(pi/180))*(spaco2+disto+disto_lab)-0.5
  y_textx1<-0.5*(sqrt(3)/2)+cos(60*(pi/180))*(spaco2+disto+disto_lab)
  x_textx2<-0.5-0.5
  y_textx2<-(-(spaco2+disto+disto_lab))
  x_textx3<-0.75+sin(60*(pi/180))*(spaco2+disto+disto_lab)-0.5
  y_textx3<-0.5*(sqrt(3)/2)+cos(60*(pi/180))*(spaco2+disto+disto_lab)
  
  #The respective position on the arrow
  arrx_x1<-0.25-sin(60*(pi/180))*(spaco2+disto)-0.5
  arry_x1<-0.5*(sqrt(3)/2)+cos(60*(pi/180))*(spaco2+disto)
  arrx_x3<-0.75+sin(60*(pi/180))*(spaco2+disto)-0.5
  arry_x3<-0.5*(sqrt(3)/2)+cos(60*(pi/180))*(spaco2+disto)
  
  #Calculate gradient triangle we need
  grad_x<-cos(60*(pi/180))*half_length_arrow
  grad_y<-sin(60*(pi/180))*half_length_arrow
  
  #Add arrows
  arrows(arrx_x1+grad_x,
         arry_x1+grad_y,
         arrx_x1-grad_x,
         arry_x1-grad_y,
         lend=3,length=0.1,lwd=arrow_lwd)
  arrows(x_textx2-half_length_arrow,
         -spaco2-disto,
         x_textx2+half_length_arrow,
         -spaco2-disto,
         lend=3,length=0.1,lwd=arrow_lwd)
  arrows(arrx_x3+grad_x,
         arry_x3-grad_y,
         arrx_x3-grad_x,
         arry_x3+grad_y,
         lend=3,length=0.1,lwd=arrow_lwd)
  
  #Add axis names
  
  safe_xpd<-par("xpd")
  par(xpd=NA)
  text(x_textx1,
       y_textx1,
       substitute(italic(b),list(b=title_names[1])),
       srt=(atan(sqrt(3)*y_over_x_ext)*180)/pi,
       cex=title_cex)
  text(x_textx2,
       y_textx2,
       substitute(italic(b),list(b=title_names[2])),srt=0,
       cex=title_cex)
  text(x_textx3,
       y_textx3,
       substitute(italic(b),list(b=title_names[3])),
       srt=360-(atan(sqrt(3)*y_over_x_ext)*180)/pi,
       cex=title_cex)
  par(xpd=safe_xpd)
  
  #CrIs
  if(length(est_coord)>0){
    
    #Go through each CrI table one by one
    for(est_list in 1:length(est_coord)){
      
      #CrI hexagons
      
      #6 intersections
      #Starting at x1 upper, x2 lower
      #x1 upper: y = -sqrt(3) * x + (1-est_coord[[est_list]][1,3])*sqrt(3)
      #x2 lower: y = sqrt(3) * x - est_coord[[est_list]][2,2]*sqrt(3)
      x_x1_upp_x2_low<-(((1-est_coord[[est_list]][1,3])*sqrt(3)+est_coord[[est_list]][2,2]*sqrt(3))/(2*sqrt(3)))-0.5
      y_x1_upp_x2_low<-sqrt(3)*(x_x1_upp_x2_low+0.5)-est_coord[[est_list]][2,2]*sqrt(3)
      
      #x1 upper, x3 lower
      #x1 upper: y = -sqrt(3) * x + (1-est_coord[[est_list]][1,3])*sqrt(3)
      #x3 lower: y = (sqrt(3)/2)*est_coord[[est_list]][3,2]
      x_x1_upp_x3_low<-(((sqrt(3)/2)*est_coord[[est_list]][3,2]-(1-est_coord[[est_list]][1,3])*sqrt(3))/(-sqrt(3)))-0.5
      y_x1_upp_x3_low<-(sqrt(3)/2)*est_coord[[est_list]][3,2]
      
      #x2 upper, x3 lower
      #x2 upper: y = sqrt(3) * x - est_coord[[est_list]][2,3]*sqrt(3)
      #x3 lower: y = (sqrt(3)/2)*est_coord[[est_list]][3,2]
      x_x2_upp_x3_low<-(((sqrt(3)/2)*est_coord[[est_list]][3,2]+est_coord[[est_list]][2,3]*sqrt(3))/(sqrt(3)))-0.5
      y_x2_upp_x3_low<-(sqrt(3)/2)*est_coord[[est_list]][3,2]
      
      #x2 upper, x1 lower
      #x2 upper: y = sqrt(3) * x - est_coord[[est_list]][2,3]*sqrt(3)
      #x1 lower: y = -sqrt(3) * x + (1-est_coord[[est_list]][1,2])*sqrt(3)
      x_x2_upp_x1_low<-(((1-est_coord[[est_list]][1,2])*sqrt(3)+est_coord[[est_list]][2,3]*sqrt(3))/(2*sqrt(3)))-0.5
      y_x2_upp_x1_low<-sqrt(3)*(x_x2_upp_x1_low+0.5)-est_coord[[est_list]][2,3]*sqrt(3)
      
      #x3 upper, x1 lower
      #x3 upper: y = (sqrt(3)/2)*est_coord[[est_list]][3,3]
      #x1 lower: y = -sqrt(3) * x + (1-est_coord[[est_list]][1,2])*sqrt(3)
      x_x3_upp_x1_low<-(((sqrt(3)/2)*est_coord[[est_list]][3,3]-(1-est_coord[[est_list]][1,2])*sqrt(3))/(-sqrt(3)))-0.5
      y_x3_upp_x1_low<-(sqrt(3)/2)*est_coord[[est_list]][3,3]
      
      #x3 upper, x2 lower
      #x3 upper: y = (sqrt(3)/2)*est_coord[[est_list]][3,3]
      #x2 lower: y = sqrt(3) * x - est_coord[[est_list]][2,2]*sqrt(3)
      x_x3_upp_x2_low<-(((sqrt(3)/2)*est_coord[[est_list]][3,3]+est_coord[[est_list]][2,2]*sqrt(3))/(sqrt(3)))-0.5
      y_x3_upp_x2_low<-(sqrt(3)/2)*est_coord[[est_list]][3,3]
      
      #Add CrI hexagon (with little black outline)
      polygon(c(x_x1_upp_x2_low,x_x1_upp_x3_low,x_x2_upp_x3_low,
                x_x2_upp_x1_low,x_x3_upp_x1_low,x_x3_upp_x2_low),
              c(y_x1_upp_x2_low,y_x1_upp_x3_low,y_x2_upp_x3_low,
                y_x2_upp_x1_low,y_x3_upp_x1_low,y_x3_upp_x2_low),
              col=adjustcolor("white",0),border="black",lwd=3.5)
      polygon(c(x_x1_upp_x2_low,x_x1_upp_x3_low,x_x2_upp_x3_low,
                x_x2_upp_x1_low,x_x3_upp_x1_low,x_x3_upp_x2_low),
              c(y_x1_upp_x2_low,y_x1_upp_x3_low,y_x2_upp_x3_low,
                y_x2_upp_x1_low,y_x3_upp_x1_low,y_x3_upp_x2_low),
              col=adjustcolor("white",0),border=est_cols[est_list],lwd=2)
    }
    
    #Go through each CrI table one by one
    for(est_list in 1:length(est_coord)){
      
      #Translate coordinates into ternary space
      est_x<-transf_x(left=est_coord[[est_list]][1,1],bottom=est_coord[[est_list]][2,1],right=est_coord[[est_list]][3,1])
      est_y<-transf_y(left=est_coord[[est_list]][1,1],bottom=est_coord[[est_list]][2,1],right=est_coord[[est_list]][3,1])
      
      #Helpful estimator lines (going from estimator to outer triangle)
      if(est_line){
        x1<-((sqrt(3)*(1-est_coord[[est_list]][1,1])-2*spaco2)/(2*sqrt(3)))-0.5
        y1<-sqrt(3)*(x1+0.5)+2*spaco2
        
        x2<-est_coord[[est_list]][2,1]-spaco2*(3^(-(1/2)))-0.5
        y2<-(-spaco2)
        
        x3<-1-(est_coord[[est_list]][3,1]/2)+(spaco2/sin(60*(pi/180)))-0.5
        y3<-est_coord[[est_list]][3,1]*(sqrt(3)/2)
        #Add lines (with little black outline)
        segments(rep(est_x,3),rep(est_y,3),
                 c(x1,x2,x3),c(y1,y2,y3),col=adjustcolor("black",0.5),lwd=2)
        segments(rep(est_x,3),rep(est_y,3),
                 c(x1,x2,x3),c(y1,y2,y3),col=est_cols[est_list],lwd=1.6)
      }
      
      #Add estimator
      points(est_x,est_y,pch=22,
             bg=est_cols[est_list],cex=1.3)
    }
    
  }
  
  #Close png if wished
  if(output_file!=""){
    invisible(dev.off())
  }
  
  #Restore old mar setting
  par(mar=old_mar)
}
```


## Figure 1 

6 Plots in one layout: 

* Blank explanatory plot
* *cydno* males
* F1 males
* *melpomene* males
* backcross to *cydno* males by genotype on chromosome 18 QTL
* backcross to *melpomene* males

We'll save this plot directly to a png. Check comments within the code chunk.

```{r}
png("Figure1.png",height=3000,width=3600,res=300)

# Set layout
layout(matrix(c(rep(1,5),rep(2,2),rep(3,5),
                rep(4,4),rep(5,4),rep(6,4),
                rep(4,4),rep(5,4),rep(6,4),
                rep(7,6),rep(8,6),
                rep(7,6),rep(8,6)),nrow=5,byrow=T))

# Plot empty plot
par(mar=c(0,0,0,0))
plot(1,1,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")

# Save old xpd settings
old_xpd<-par("xpd")
# Set xpd to FALSE
par(xpd=FALSE)

# Plot explanatory graph
ternary_choice(
  output_file="",
  data_fr=c(),
  est_coord=c(),
  cont_yn=F,
  ticks_inner_outer="outer",
  spaco2=0.07,
  disto=0.1,
  half_length_arrow=0.3,
  arrow_lwd=2,
  disto_lab=0.1,
  tick_length1=0,
  space_around=c(((1-((sqrt(3))/2)))+0.14,0.25,0.36,0.25),
  line_at=seq(0,1,0.1),
  tck_at=seq(0,1,0.1),
  axis_cex=0,
  title_names=c("Only cydno","Only melpomene","Both"),
  title_cex=1.8,
  scales_to_plot=2
)

# Add photos of butterflies
par(xpd=NA)
# cydno
# height to width ratio of cydno photo
cyd_rat<-dim(cydno)[1]/dim(cydno)[2]
# add cydno
rasterImage(cydno,-0.97,-0.15,-0.64,-0.15+0.33*cyd_rat)
# melpomene
# height to width ratio of melpomene photo
mel_rat<-dim(melpomene)[1]/dim(melpomene)[2]
# add melpomene
rasterImage(melpomene,0.64,-0.15,0.97,-0.15+0.33*mel_rat)

# Both
rasterImage(cydno,-0.34,1.01,-0.01,1.01+0.33*cyd_rat)
rasterImage(melpomene,0.01,1.01,0.34,1.01+0.33*mel_rat)

# Inform about cydno-melpomene scale
text(-0.21,0.15,
     substitute(paste(italic('cyd.'),"+")),cex=1.7)
text(0.21,0.16,
     substitute(paste(italic('mel.'),"+")),cex=1.7)

# Put dot sizes
text(1,1.05,"# Trials",cex=1.5)
points(rep(0.93,5),seq(0.2,(sqrt(3))/2,length.out=5),pch=21,bg=adjustcolor("black",0.2),cex=1.5*sqrt(1:5))
text(rep(1.08,5),seq(0.2,(sqrt(3))/2,length.out=5),as.character(1:5),cex=1.5)

par(xpd=F)

# Plot empty plot
plot(1,1,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")


# Plot cydno
ternary_choice(
  output_file="",
  data_fr=tern_pref[tern_pref$Type=="CP",c("prop_trials_cydno_court_only","prop_trials_melp_court_only","prop_trials_court_both")],
  dot_cols=rep("dodgerblue3",sum(tern_pref$Type=="CP")),
  dot_sizes=1.5*sqrt(tern_pref$total_trials_with_response[tern_pref$Type=="CP"]),
  transp=0.2,
  jitter_dots=3,
  est_coord=estimator_table_1[1],
  est_cols="dodgerblue3",
  est_line=T,
  cont_yn=F,
  remove_outer=F,
  ticks_inner_outer="outer",
  spaco2=0.07,
  disto=100,
  tick_length1=0.02,
  space_around=c(((1-((sqrt(3))/2)))+0.08,0.16,0.24,0.16),
  line_at=seq(0,1,0.1),
  tck_at=seq(0,1,0.1),
  axis_l=seq(0,1,0.1),
  axis_cex=1.2,
  axis_d1=0.04,
  title_names=c("","",""),
  scales_to_plot=2
)

# Legend:
par(xpd=NA)
text(-0.35,0.85,substitute(italic(a),list(a="cydno")),cex=2,pos=2,offset=0)
par(xpd=F)


# Plot F1
ternary_choice(
  output_file="",
  data_fr=tern_pref[tern_pref$Type=="CPxMP",c("prop_trials_cydno_court_only","prop_trials_melp_court_only","prop_trials_court_both")],
  dot_cols=rep("orange",sum(tern_pref$Type=="CPxMP")),
  dot_sizes=1.5*sqrt(tern_pref$total_trials_with_response[tern_pref$Type=="CPxMP"]),
  jitter_dots=3,
  transp=0.2,
  est_coord=estimator_table_1[3],
  est_cols="orange",
  est_line=T,
  cont_yn=F,
  remove_outer=F,
  ticks_inner_outer="outer",
  spaco2=0.07,
  disto=100,
  tick_length1=0.02,
  space_around=c(((1-((sqrt(3))/2)))+0.08,0.16,0.24,0.16),
  line_at=seq(0,1,0.1),
  tck_at=seq(0,1,0.1),
  axis_l=seq(0,1,0.1),
  axis_cex=1.2,
  axis_d1=0.04,
  title_names=c("","",""),
  scales_to_plot=2
)

# Legend:
par(xpd=NA)
text(-0.35,0.85,"F1",cex=2,pos=2,offset=0)
par(xpd=F)


# Plot melpomene
ternary_choice(
  output_file="",
  data_fr=tern_pref[tern_pref$Type=="MP",c("prop_trials_cydno_court_only","prop_trials_melp_court_only","prop_trials_court_both")],
  dot_cols=rep("orange",sum(tern_pref$Type=="MP")),
  dot_sizes=1.5*sqrt(tern_pref$total_trials_with_response[tern_pref$Type=="MP"]),
  jitter_dots=3,
  transp=0.2,
  est_coord=estimator_table_1[2],
  est_cols="orange",
  est_line=T,
  cont_yn=F,
  remove_outer=F,
  ticks_inner_outer="outer",
  spaco2=0.07,
  disto=100,
  tick_length1=0.02,
  space_around=c(((1-((sqrt(3))/2)))+0.08,0.16,0.24,0.16),
  line_at=seq(0,1,0.1),
  tck_at=seq(0,1,0.1),
  axis_l=seq(0,1,0.1),
  axis_cex=1.2,
  axis_d1=0.04,
  title_names=c("","",""),
  scales_to_plot=2
)

# Legend:
par(xpd=NA)
text(-0.35,0.85,substitute(italic(a),list(a="melpomene")),cex=2,pos=2,offset=0)
par(xpd=F)


# Plot cydno BC by genotype on chromosome 18 peak
ternary_choice(
  output_file="",
  data_fr=tern_pref[tern_pref$Type=="CPx(CPxMP)",c("prop_trials_cydno_court_only","prop_trials_melp_court_only","prop_trials_court_both")],
  dot_cols=sapply((1:nrow(tern_pref_cydBC))[tern_pref_cydBC$geno18%in%c("A","B")],function(x) ifelse(tern_pref_cydBC$geno18[x]=="A","dodgerblue3","orange")),
  dot_sizes=1.5*sqrt(tern_pref$total_trials_with_response[tern_pref$Type=="CPx(CPxMP)"]),
  jitter_dots=3,
  transp=0.2,
  est_coord=estimator_table_2_geno18,
  est_cols=c("dodgerblue3","orange"),
  est_line=T,
  cont_yn=F,
  remove_outer=F,
  ticks_inner_outer="outer",
  spaco2=0.07,
  disto=100,
  tick_length1=0.02,
  space_around=c(((1-((sqrt(3))/2)))+0.08,0.82,0.24,0.16), #Give more space on the left side to slide the triangle towards center
  line_at=seq(0,1,0.1),
  tck_at=seq(0,1,0.1),
  axis_l=seq(0,1,0.1),
  axis_cex=1.2,
  axis_d1=0.04,
  title_names=c("","",""),
  scales_to_plot=2
)

# Legend:
par(xpd=NA)
text(-0.35,0.75,substitute(paste("BC to ",italic('cyd.'))),cex=2,pos=2,offset=0)
par(xpd=F)


# Plot melpomene backcross
ternary_choice(
  output_file="",
  data_fr=tern_pref[tern_pref$Type=="MPx(CPxMP)",c("prop_trials_cydno_court_only","prop_trials_melp_court_only","prop_trials_court_both")],
  dot_cols=rep("orange",sum(tern_pref$Type=="MPx(CPxMP)")),
  dot_sizes=1.5*sqrt(tern_pref$total_trials_with_response[tern_pref$Type=="MPx(CPxMP)"]),
  jitter_dots=3,
  transp=0.2,
  est_coord=estimator_table_1[5],
  est_cols="orange",
  est_line=T,
  cont_yn=F,
  remove_outer=F,
  ticks_inner_outer="outer",
  spaco2=0.07,
  disto=100,
  tick_length1=0.02,
  space_around=c(((1-((sqrt(3))/2)))+0.08,0.16,0.24,0.82), #Give more space on the right side to slide the triangle towards center
  line_at=seq(0,1,0.1),
  tck_at=seq(0,1,0.1),
  axis_l=seq(0,1,0.1),
  axis_cex=1.2,
  axis_d1=0.04,
  title_names=c("","",""),
  scales_to_plot=2
)

# Legend:
par(xpd=NA)
text(-0.35,0.75,substitute(paste("BC to ",italic('mel.'))),cex=2,pos=2,offset=0)
par(xpd=F)

invisible(dev.off())
```

Look at png we just made
```{r,fig.width=8.4, fig.height=7}
fig1<-readPNG("Figure1.png")
grid.raster(fig1)
```

For completeness, we look at the estimators for backcross-to-*cydno* males independent of their genotype:

```{r}
cyd_BC<-estimator_table_1[[4]]
row.names(cyd_BC)<-c("cyd_only","melp_only","both")
print(cyd_BC)
```

## Supplementary Figure 1

5 Plots in one layout: 

* Blank explanatory plot
* backcross to *cydno* males homozygous on chromosome 1 & 18 QTL peak
* backcross to *cydno* males homozygous on chromosome 1 & heterzygous on chromosome 18 QTL peak
* backcross to *cydno* males heterozygous on chromosome 1 & homozygous on chromosome 18 QTL peak
* backcross to *cydno* males heterozygous on chromosome 1 & 18 QTL peak

Safe again as a png

```{r}
# Resize png according to omi settings
png("Suppl_Figure1.png",height=(10/9.5)*3000,width=(8/7.5)*2400,res=300)

# Set layout
layout(matrix(c(rep(1,9),rep(2,6),rep(3,9),
                rep(4,12),rep(5,12),
                rep(4,12),rep(5,12),
                rep(6,12),rep(7,12),
                rep(6,12),rep(7,12)),nrow=5,byrow=T))

#Give outer margins for later adding the labels.
par(omi=c(0.5,0.5,0,0))

# Plot empty plot
par(mar=c(0,0,0,0))
plot(1,1,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")

# Plot explanatory graph
ternary_choice(
  output_file="",
  data_fr=c(),
  est_coord=c(),
  cont_yn=F,
  ticks_inner_outer="outer",
  spaco2=0.07,
  disto=0.1,
  half_length_arrow=0.3,
  arrow_lwd=2,
  disto_lab=0.1,
  tick_length1=0,
  space_around=c(((1-((sqrt(3))/2)))+0.14,0.25,0.36,0.25),
  line_at=seq(0,1,0.1),
  tck_at=seq(0,1,0.1),
  axis_cex=0,
  title_names=c("Only cydno","Only melpomene","Both"),
  title_cex=1.8,
  scales_to_plot=2
)

# Add photos of butterflies
par(xpd=NA)
# cydno
# height to width ratio of cydno photo
cyd_rat<-dim(cydno)[1]/dim(cydno)[2]
# add cydno
rasterImage(cydno,-0.97,-0.15,-0.64,-0.15+0.33*cyd_rat)
# melpomene
# height to width ratio of melpomene photo
mel_rat<-dim(melpomene)[1]/dim(melpomene)[2]
# add melpomene
rasterImage(melpomene,0.64,-0.15,0.97,-0.15+0.33*mel_rat)

# Both
rasterImage(cydno,-0.34,1.01,-0.01,1.01+0.33*cyd_rat)
rasterImage(melpomene,0.01,1.01,0.34,1.01+0.33*mel_rat)

# Inform about scale
text(-0.21,0.15,
     substitute(paste(italic('cyd.'),"+")),cex=1.7)
text(0.21,0.16,
     substitute(paste(italic('mel.'),"+")),cex=1.7)

# Put dot sizes
text(1,1.05,"# Trials",cex=1.5)
points(rep(0.93,5),seq(0.2,(sqrt(3))/2,length.out=5),pch=21,bg=adjustcolor("black",0.2),cex=1.5*sqrt(1:5))
text(rep(1.08,5),seq(0.2,(sqrt(3))/2,length.out=5),as.character(1:5),cex=1.5)

par(xpd=F)

# Plot empty plot
plot(1,1,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")


# Now all cydno backcross males homozygous on geno1 and geno18
ternary_choice(
  output_file="",
  data_fr=tern_pref_cydBC[tern_pref_cydBC$geno1=="A"&tern_pref_cydBC$geno18=="A",c("prop_trials_cydno_court_only","prop_trials_melp_court_only","prop_trials_court_both")],
  dot_cols=rep("black",sum(tern_pref_cydBC$geno1=="A"&tern_pref_cydBC$geno18=="A")),
  dot_sizes=1.5*sqrt(tern_pref_cydBC$total_trials_with_response[tern_pref_cydBC$geno1=="A"&tern_pref_cydBC$geno18=="A"]),
  jitter_dots=3,
  transp=0.2,
  est_coord=estimator_table_2_geno18_with_geno1_homo[1],
  est_cols=c("black"),
  est_line=T,
  cont_yn=F,
  remove_outer=F,
  ticks_inner_outer="outer",
  spaco2=0.07,
  disto=100,
  tick_length1=0.02,
  space_around=c(((1-((sqrt(3))/2)))+0.08,0.16,0.24,0.16),
  line_at=seq(0,1,0.1),
  tck_at=seq(0,1,0.1),
  axis_l=seq(0,1,0.1),
  axis_cex=1,
  axis_d1=0.035,
  title_names=c("","",""),
  scales_to_plot=c(2)
)

par(xpd=NA)
text(-0.75,0.5,expression('chr1'^ italic('cyd:cyd')),cex=2.3,offset=0,srt=90)
par(xpd=F)

# Now all cydno backcross males homozygous on geno1 and heterozygous on geno18
ternary_choice(
  output_file="",
  data_fr=tern_pref_cydBC[tern_pref_cydBC$geno1=="A"&tern_pref_cydBC$geno18=="B",c("prop_trials_cydno_court_only","prop_trials_melp_court_only","prop_trials_court_both")],
  dot_cols=rep("black",sum(tern_pref_cydBC$geno1=="A"&tern_pref_cydBC$geno18=="B")),
  dot_sizes=1.5*sqrt(tern_pref_cydBC$total_trials_with_response[tern_pref_cydBC$geno1=="A"&tern_pref_cydBC$geno18=="B"]),
  jitter_dots=3,
  transp=0.2,
  est_coord=estimator_table_2_geno18_with_geno1_homo[2],
  est_cols=c("black"),
  est_line=T,
  cont_yn=F,
  remove_outer=F,
  ticks_inner_outer="outer",
  spaco2=0.07,
  disto=100,
  tick_length1=0.02,
  space_around=c(((1-((sqrt(3))/2)))+0.08,0.16,0.24,0.16),
  line_at=seq(0,1,0.1),
  tck_at=seq(0,1,0.1),
  axis_l=seq(0,1,0.1),
  axis_cex=1,
  axis_d1=0.035,
  title_names=c("","",""),
  scales_to_plot=c(2)
)

# Now all cydno backcross males heterozygous on geno1 and homozygous on geno18
ternary_choice(
  output_file="",
  data_fr=tern_pref_cydBC[tern_pref_cydBC$geno1=="B"&tern_pref_cydBC$geno18=="A",c("prop_trials_cydno_court_only","prop_trials_melp_court_only","prop_trials_court_both")],
  dot_cols=rep("black",sum(tern_pref_cydBC$geno1=="B"&tern_pref_cydBC$geno18=="A")),
  dot_sizes=1.5*sqrt(tern_pref_cydBC$total_trials_with_response[tern_pref_cydBC$geno1=="B"&tern_pref_cydBC$geno18=="A"]),
  jitter_dots=3,
  transp=0.2,
  est_coord=estimator_table_2_geno18_with_geno1_hetero[1],
  est_cols=c("black"),
  est_line=T,
  cont_yn=F,
  remove_outer=F,
  ticks_inner_outer="outer",
  spaco2=0.07,
  disto=100,
  tick_length1=0.02,
  space_around=c(((1-((sqrt(3))/2)))+0.08,0.16,0.24,0.16),
  line_at=seq(0,1,0.1),
  tck_at=seq(0,1,0.1),
  axis_l=seq(0,1,0.1),
  axis_cex=1,
  axis_d1=0.035,
  title_names=c("","",""),
  scales_to_plot=c(2)
)

par(xpd=NA)
text(-0.75,0.5,expression('chr1'^ italic('cyd:mel')),cex=2.3,offset=0,srt=90)
par(xpd=F)

par(xpd=NA)
text(0,-0.3,expression('chr18'^ italic('cyd:cyd')),cex=2.3,offset=0)
par(xpd=F)

# Now all cydno backcross males heterozygous on geno1 and geno18
ternary_choice(
  output_file="",
  data_fr=tern_pref_cydBC[tern_pref_cydBC$geno1=="B"&tern_pref_cydBC$geno18=="B",c("prop_trials_cydno_court_only","prop_trials_melp_court_only","prop_trials_court_both")],
  dot_cols=rep("black",sum(tern_pref_cydBC$geno1=="B"&tern_pref_cydBC$geno18=="B")),
  dot_sizes=1.5*sqrt(tern_pref_cydBC$total_trials_with_response[tern_pref_cydBC$geno1=="B"&tern_pref_cydBC$geno18=="B"]),
  jitter_dots=3,
  transp=0.2,
  est_coord=estimator_table_2_geno18_with_geno1_hetero[2],
  est_cols=c("black"),
  est_line=T,
  cont_yn=F,
  remove_outer=F,
  ticks_inner_outer="outer",
  spaco2=0.07,
  disto=100,
  tick_length1=0.02,
  space_around=c(((1-((sqrt(3))/2)))+0.08,0.16,0.24,0.16),
  line_at=seq(0,1,0.1),
  tck_at=seq(0,1,0.1),
  axis_l=seq(0,1,0.1),
  axis_cex=1,
  axis_d1=0.035,
  title_names=c("","",""),
  scales_to_plot=c(2)
)

par(xpd=NA)
text(0,-0.3,expression('chr18'^ italic('cyd:mel')),cex=2.3,offset=0)
par(xpd=F)

invisible(dev.off())

#Restore old xpd settings
par(xpd=old_xpd)
```

Look at png we just made
```{r,fig.width=6, fig.height=7.4}
suppl_fig1<-readPNG("Suppl_Figure1.png")
grid.raster(suppl_fig1)
```
