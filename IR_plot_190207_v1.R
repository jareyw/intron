library(reshape2)
library(ggplot2)
library(gdata)
library(plyr)
library(readr)

get_bg_colors = function(bg_color){
  if(bg_color=="white"){
    front_color="black"
    back_color="white"
  }
  if(bg_color=="black"){
    front_color="white"
    back_color="black"
  }
  return(c(front_color,back_color))
}

## SOURCE: https://stackoverflow.com/questions/28273716/r-implementation-for-finding-the-longest-common-starting-substrings-in-a-set-of
comsub<-function(x) {
  # sort the vector
  x<-sort(x)
  # split the first and last element by character
  d_x<-strsplit(x[c(1,length(x))],"")
  # search for the first not common element and so, get the last matching one
  der_com<-match(FALSE,do.call("==",d_x))-1
  # if there is no matching element, return an empty vector, else return the common part
  ifelse(der_com==0,return(character(0)),return(substr(x[1],1,der_com)))
}

###
###ONLY REPLETE FUNCTION
###
###RE-WRITE TO BE MORE INTUITIVE (LESS HARD CODING)
###

get_IR_table = function(data_table,sample_names,color_vector,average){
  #Isolate from data table only the specified samples
  new_colnames = grep(paste(sample_names,collapse="|"),colnames(data_table))
  data_table = cbind(data_table[,1:2],data_table[,new_colnames])
  k = colnames(data_table)[-c(1,2)]
  
  print(new_colnames)
  print(k)
  
  if(average==FALSE){
    k_labs=k
    if(length(color_vector)==2){
      color_vector = colorRampPalette(color_vector)(length(k))
    }
    else(color_vector = color_vector)
  }
  if(average==TRUE){
    new_data_table = data_table[,1:2]
    k_labs = sample_names
    
    for(sample in sample_names){
      temp_table = data_table[,grep(sample,colnames(data_table))]
      new_data_table[[sample]] = rowMeans(temp_table)
    }
    n=1
    data_table = new_data_table
  }
  print(color_vector)
  
  output = list(data_table,k_labs,color_vector)
  return(output)
}

exact_ks_test = function(data_file,sample_names,color_vector,bg_color,n,average=FALSE){
  data_table = read_delim(data_file,"\t",escape_double = FALSE,trim_ws = TRUE)
  front_color = get_bg_colors(bg_color)[1]
  back_color = get_bg_colors(bg_color)[2]
  data_table = as.data.frame(data_table)
  
  data_out = get_IR_table(data_table,sample_names,color_vector,average)
  data_table = as.data.frame(data_out[[1]])
  k_labs = data_out[[2]]
  color_vector = data_out[[3]]
  
  #need to remove rows (i.e. for all conditions) in which there are 0-values if the plot is in log scale
  print(nrow(data_table))
  for (i in k_labs){
    data_table=data_table[with(data_table,eval(as.name(i)))>0,]
  }
  #data_table=data_table[with(data_table,eval(as.name(x))>0 & eval(as.name(y))>0),]
  print(nrow(data_table)) #spot check to see how many rows you are losing
  
  sample_values_1 = data_table[[sample_names[1]]]
  sample_values_2 = data_table[[sample_names[2]]]
  print(ks.test(sample_values_1,sample_values_2,alternative="greater",exact=TRUE)$p.value)
}

plot_logIR_eCDF = function(data_file,sample_names,color_vector,bg_color,n,
                           x_min,x_max=0,output,average=FALSE,pseudo=FALSE,legend_pos=FALSE){
  data_table = read_delim(data_file,"\t",escape_double = FALSE,trim_ws = TRUE)
  front_color = get_bg_colors(bg_color)[1]
  back_color = get_bg_colors(bg_color)[2]
  data_table = as.data.frame(data_table)
  
  data_out = get_IR_table(data_table,sample_names,color_vector,average)
  data_table = as.data.frame(data_out[[1]])
  k_labs = data_out[[2]]
  color_vector = data_out[[3]]
  
  if(pseudo!=FALSE){
    data_table[,3:4] = data_table[,3:4]+pseudo
  }

  #need to remove rows (i.e. for all conditions) in which there are 0-values if the plot is in log scale
  print(nrow(data_table))
  for (i in k_labs){
    data_table=data_table[with(data_table,eval(as.name(i)))>0,]
  }
  #data_table=data_table[with(data_table,eval(as.name(x))>0 & eval(as.name(y))>0),]
  print(nrow(data_table)) #spot check to see how many rows you are losing
  
  sample_values_1 = data_table[[sample_names[1]]]
  sample_values_2 = data_table[[sample_names[2]]]
  print(ks.test(sample_values_1,sample_values_2,alternative="greater",exact=TRUE)$p.value)
  print(ks.test(sample_values_1,sample_values_2,alternative="less",exact=TRUE)$p.value)
  
  long_data = melt(data_table,id.vars=c("intron_id","info"),variable.name="Condition",value.name="IR_Ratio")
  sub_data = long_data[long_data$Condition %in% k_labs,]
  print(nrow(sub_data))
  
  sub_data$log_IR_Ratio = log10(sub_data$IR_Ratio)
  print(k_labs)
  
  if(legend_pos!=FALSE){
    legend_pos = as.vector(legend_pos)
  }else{legend_pos="none"}
  
  ggplot(sub_data,aes(x=log_IR_Ratio,color=Condition))+
    stat_ecdf(size=1.2)+
    ggtitle(paste(nrow(data_table),"introns",sep=" "))+
    ylab("Fraction of introns")+
    xlab("Log10(IR Score)")+
    coord_cartesian(xlim=c(x_min,x_max))+
    scale_color_manual(values=color_vector,labels=k_labs)+
    theme(panel.border=element_blank(),
          panel.background=element_blank(),
          plot.background=element_rect(fill="transparent",colour=NA),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.ticks=element_line(colour=front_color),
          axis.title=element_text(colour=front_color,size=18),
          axis.text=element_text(colour=front_color,size=18),
          legend.background = element_rect(fill="transparent"),
          legend.title = element_text(color=front_color,size=16),
          legend.key=element_rect(fill="transparent",color=NA),
          legend.text=element_text(color=front_color,size=12),
          legend.justification=c(1,0),
          legend.position=legend_pos,
          plot.title=element_text(color=front_color))+
    theme(axis.line.x=element_line(color=front_color))+
    theme(axis.line.y=element_line(color=front_color))+
    theme(panel.border = element_rect(color=front_color,fill=NA,size=2))
  
  #Most images in the past have been width = 6, height = 6
  ggsave(paste(output,k_labs[1],"_",k_labs[2],"_",bg_color,"_eCDF.pdf",sep=""),width=6,height=6,bg="transparent")
  ggsave(paste(output,k_labs[1],"_",k_labs[2],"_",bg_color,"_eCDF.png",sep=""),width=6,height=6,bg="transparent")
  dev.off()
}

plot_logIR_duoScatter = function(data_file,sample_names,color_vector,bg_color,n,x_min,x_max=0,output,title="",average=FALSE,pseudo=FALSE){
  data_table = read_delim(data_file,"\t",escape_double = FALSE,trim_ws = TRUE)
  front_color = get_bg_colors(bg_color)[1]
  back_color = get_bg_colors(bg_color)[2]
  data_table = as.data.frame(data_table)
  
  data_out = get_IR_table(data_table,sample_names,color_vector,average)
  data_table = as.data.frame(data_out[[1]])
  k_labs = data_out[[2]]
  color_vector = data_out[[3]]
  
  if(pseudo!=FALSE){
    data_table[,3:4] = data_table[,3:4]+pseudo
  }
  
  #need to remove rows (i.e. for all conditions) in which there are 0-values if the plot is in log scale
  print(nrow(data_table))
  for (i in k_labs){
    data_table=data_table[with(data_table,eval(as.name(i)))>0,]
  }
  #data_table=data_table[with(data_table,eval(as.name(x))>0 & eval(as.name(y))>0),]
  print(nrow(data_table)) #spot check to see how many rows you are losing
  #data_table[,3:4] = log10(data_table[,3:4])
  
  #If you did the stupid thing and named your factors starting with non-characters
  # k_labs[1] = paste("`",k_labs[1],"`",sep="",collapse="")
  # k_labs[2] = paste("`",k_labs[2],"`",sep="",collapse="")
  
  ggplot(data_table,aes_string(x=k_labs[2],y=k_labs[1]))+
    stat_density2d(aes(fill=..level..,alpha=..level..),
                   geom='polygon',
                   color='grey50')+
    scale_fill_gradient(low="white",high="red")+
    scale_alpha(range=c(0,1),guide=F)+
    geom_point(color=front_color,alpha=0.1)+
    ggtitle(paste(nrow(data_table),"introns",sep=" "))+
    xlab(paste("Log10(IR Score) ",k_labs[2],sep=""))+
    ylab(paste("Log10(IR Score) ",k_labs[1],sep=""))+
    xlim(x_min,x_max)+
    ylim(x_min,x_max)+
    geom_abline(intercept=0,slope=1,color="blue",linetype="dashed",size=1.5)+
    theme(panel.border=element_blank(),
          panel.background=element_blank(),
          plot.background=element_rect(fill="transparent",colour=NA),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.ticks=element_line(colour=front_color),
          axis.title=element_text(colour=front_color,size=18),
          axis.text=element_text(colour=front_color,size=18),
          legend.background = element_rect(fill="transparent"),
          legend.title = element_text(color=front_color,size=16),
          legend.key=element_rect(fill="transparent",color=NA),
          legend.text=element_text(color=front_color,size=12),
          legend.justification=c(1,0),
          legend.position=c(0.3,0.6),
          plot.title=element_text(color=front_color))+
    theme(axis.line.x=element_line(color=front_color))+
    theme(axis.line.y=element_line(color=front_color))+
    theme(panel.border = element_rect(color=front_color,fill=NA,size=2))
  
  #Most images in the past have been width = 6, height = 6
  ggsave(paste(output,k_labs[1],"_",k_labs[2],"_",bg_color,"_",title,"_scatter.png",sep=""),width=6,height=6,bg="transparent")
  dev.off()
}

plot_logIR_FC_density = function(df,output,title,bg_color){
  front_color = get_bg_colors(bg_color)[1]
  back_color = get_bg_colors(bg_color)[2]
  
  max = 2.5
  df$DMSO_FC[df$DMSO_FC > 2.5] <- 2.5
  df$`8800_FC`[df$`8800_FC` > 2.5] <- 2.5
  # df$INPUT_FC[df$INPUT_FC > 10^2.5] <- 10^2.5
  # df$J2_FC[df$J2_FC > 10^2.5] <- 10^2.5
  
  ggplot(data=df,aes(log10(INPUT_FC),log10(J2_FC)))+
    #stat_density2d(aes(fill=..level..,alpha=..level..),
    #               geom='polygon',
    #               color='grey50')+
    #scale_fill_gradient(low="white",high="red")+
    #scale_alpha(range=c(0,1),guide=F)+
    # geom_point(alpha=0.1,size=0.3)+
    # geom_point(data = subset(df,INPUT_FC/J2_FC > 2), size = 1.25, alpha = 0.75, col="blue")+
    # geom_point(data = subset(df,INPUT_FC/J2_FC < 0.5), size = 1.25, alpha = 0.75, col="red")+
    # geom_abline(slope=1,intercept=0,color="black",linetype=2,size=1)+
    geom_point(data = subset(df, `8800_FC` > 1 & `8800_FC`>DMSO_FC), size=1, alpha=0.6, col="red")+
    geom_point(data = subset(df, DMSO_FC > 1 & DMSO_FC>`8800_FC`), size=1, alpha=0.6, col="blue")+
    geom_hline(yintercept=1,color="grey",linetype=2,size=0.5)+
    geom_vline(xintercept=1,color="grey",linetype=2,size=0.5)+
    ylim(-1,max)+
    xlim(-1,max)+
    theme_classic()+
    theme(panel.border = element_rect(color=front_color,fill=NA,size=1),
          axis.ticks=element_line(colour=front_color),
          axis.text=element_text(colour=front_color,size=16),
          axis.title=element_text(colour=front_color,size=16))
  ggsave(paste(output,title,"_","_",bg_color,"_heatscatter.png",sep=""),width=5,height=5,bg="transparent")
  dev.off()
}

log_IR_boxplots = function(data_file,sample_names,color_vector,bg_color,n,output,title="",y_label="",log=FALSE,average=FALSE,pseudo=FALSE,subset=FALSE){
  data_table = read_delim(data_file,"\t",escape_double = FALSE,trim_ws = TRUE)
  front_color = get_bg_colors(bg_color)[1]
  back_color = get_bg_colors(bg_color)[2]
  data_table = as.data.frame(data_table)
  
  data_out = get_IR_table(data_table,sample_names,color_vector,average)
  data_table = as.data.frame(data_out[[1]])
  k_labs = data_out[[2]]
  color_vector = data_out[[3]]
  
  if(pseudo!=FALSE){
    data_table[,3:4] = data_table[,3:4]+pseudo
  }
  if(subset!=FALSE){
    subset = as.vector(subset)
    data_table = data_table[subset,]
  }
  
  print(data_table)
  
  #need to remove rows (i.e. for all conditions) in which there are 0-values if the plot is in log scale
  print(nrow(data_table))
  for (i in k_labs){
    data_table=data_table[with(data_table,eval(as.name(i)))>0,]
  }
  #data_table=data_table[with(data_table,eval(as.name(x))>0 & eval(as.name(y))>0),]
  print(nrow(data_table)) #spot check to see how many rows you are losing
  
  long_data = melt(data_table,id.vars=c("intron_id","info"),variable.name="Condition",value.name="IR_Ratio")
  sub_data = long_data[long_data$Condition %in% k_labs,]
  print(nrow(sub_data))
  
  if(log!=FALSE){
    sub_data$IR_Ratio = log10(sub_data$IR_Ratio)
  }

  print(k_labs)
  
  ggplot(sub_data,aes(y=IR_Ratio,x=Condition,color=Condition))+
    geom_boxplot()+
    ggtitle(paste(nrow(data_table),"introns",sep=" "))+
    ylab(y_label)+
    scale_color_manual(values=color_vector,labels=k_labs)+
    theme(panel.border=element_blank(),
        panel.background=element_blank(),
        plot.background=element_rect(fill="transparent",colour=NA),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks=element_line(colour=front_color),
        axis.title=element_text(colour=front_color,size=18),
        axis.text=element_text(colour=front_color,size=18),
        legend.background = element_rect(fill="transparent"),
        legend.title = element_text(color=front_color,size=16),
        legend.key=element_rect(fill="transparent",color=NA),
        legend.text=element_text(color=front_color,size=12),
        legend.justification=c(1,0),
        legend.position="none",
        plot.title=element_text(color=front_color))+
    theme(axis.line.x=element_line(color=front_color))+
    theme(axis.line.y=element_line(color=front_color))+
    theme(panel.border = element_rect(color=front_color,fill=NA,size=2))
  
  #Most images in the past have been width = 6, height = 6
  ggsave(paste(output,k_labs[1],"_",k_labs[2],"_",bg_color,"_",title,"_boxplots.png",sep=""),width=6,height=6,bg="transparent")
  dev.off()
}
