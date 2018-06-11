##plot 
rm(list=ls())
setwd("E:/Paper/TP-SSU-SKAT")
library(dplyr)
library(tidyr)
library(ggplot2)
#type-I-error
library(xlsx)
dat_t = read.xlsx(file = "tpSSU_result_cov -plot.xlsx",1,header = FALSE)
type1_10 = as.matrix(dat_t[3:4,2:5])
type1_20 = as.matrix(dat_t[7:8,2:5])

type = type1_20
ty_mat = apply(type,2,as.numeric)

ty = data.frame(row.names = 1:8)
ty$method = as.factor(rep(c("tpSSU","tpSKAT","SSU","SKAT"),2))
ty$cor = as.factor(c(rep("CS",4),rep("AR-1",4)))
ty$value = c(ty_mat[1,],ty_mat[2,])

 #画布大小
                         
ggplot(ty,aes(x=factor(cor,levels = levels(ty$cor)[c(2,1)]), 
              y=value, fill=factor(method,levels = levels(ty$method)[c(4,3,2,1)])))+
  geom_bar(stat='identity',
           position=position_dodge())+ #本来不用reverse是从上到下，反过来
  scale_fill_discrete(name='legend', #图例项（或者用scale_fill_discrete)
                      labels=levels(ty$method)[c(4,3,2,1)])+
  guides(fill=guide_legend(title=NULL))+
  # scale_fill_manual(breaks = c("tpSSU", "tpSKAT", "SSU","SKAT"),
  #                   values=c("blue","red","green","violet"))+
  theme(legend.position="right")+#图例位置
  theme(legend.text = element_text(size = 10))+
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14))+
  scale_x_discrete(name='Correlation Struction', #x轴坐标名称
                   labels= levels(ty$cor)[c(2,1)])+ #离散的标签
  scale_y_continuous(name='Type I Error', #y轴坐标名称
                     breaks=seq(0,0.08,0.01),
                     limits=c(0,0.08))#连续的标签和坐标轴

# power
dat_rec = read.xlsx(file ="tpSSU_result_cov -plot.xlsx",2)
#s1
##m=10
rec_10_1c = as.matrix(dat_rec[3:10,3:7])
rec_10_1a = as.matrix(dat_rec[11:18,3:7])
##m=20
rec_20_1c = as.matrix(dat_rec[62:69,3:7])
rec_20_1a = as.matrix(dat_rec[71:78,3:7])

#s2
##m=10
rec_10_2c = as.matrix(dat_rec[22:26,3:7])
rec_10_2a = as.matrix(dat_rec[27:31,3:7])
##m=20
rec_20_2c = as.matrix(dat_rec[81:85,3:7])
rec_20_2a = as.matrix(dat_rec[86:90,3:7])
#s3
##m=10
rec_10_3c = as.matrix(dat_rec[37:44,3:7])
rec_10_3a = as.matrix(dat_rec[47:54,3:7])
##m=20
rec_20_3c = as.matrix(dat_rec[96:103,3:7])
rec_20_3a = as.matrix(dat_rec[106:113,3:7])

#s4
dat_dom = read.xlsx(file = "tpSSU_result_cov -plot.xlsx",sheetIndex = 4)
##m=10
rec_10_4c = as.matrix(dat_dom[3:7,3:7])
rec_10_4a = as.matrix(dat_dom[8:12,3:7])
##m=20
rec_20_4c = as.matrix(dat_dom[19:23,3:7])
rec_20_4a = as.matrix(dat_dom[24:28,3:7])

#s5
dat_add = read.xlsx(file="tpSSU_result_cov -plot.xlsx",3)
##m=10
rec_10_5c = as.matrix(dat_add[4:8,3:7])
rec_10_5a = as.matrix(dat_add[10:14,3:7])
##m=20
rec_20_5c = as.matrix(dat_add[20:24,3:7])
rec_20_5a = as.matrix(dat_add[26:30,3:7])


cs = rec_20_5c 
ar1=rec_20_5a
a= plotpower(cs,"CS")
b = plotpower(ar1,"AR-1")
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(a, vp = vplayout(1,1))   ###将（1,1)和(1,2)的位置画图c
print(b, vp = vplayout(1,2))   ###将(2,1)的位置画图b
# dev.off()


plotpower<-function(data,title){
  dat = as.data.frame(apply(data,2,as.numeric))
  colnames(dat) = c("beta","tpSSU","tpSKAT","SSU","SKAT")
  dat$beta = factor(round(dat$beta,2))
  dat$id = 1:nrow(dat)
  
  df1= dat %>% gather("item",value,-c("id","beta")) %>% 
    bind_cols(data.frame(item_id=rep(1:4,each=nrow(dat))))
  ##折线图
  df1$item = as.factor(df1$item)
  library(grid)
  a<-ggplot(df1,aes(id,value,
                    linetype =factor(item,levels = levels(df1$item)[c(4,3,2,1)]),
                    colour=factor(item,levels = levels(df1$item)[c(4,3,2,1)]),
                    shape =factor(item,levels = levels(df1$item)[c(4,3,2,1)]) ))+
    geom_line(size=0.5,position=position_dodge(0))+
    theme(legend.title=element_blank())+
    geom_point(size=2,position=position_dodge(0))+
    scale_x_continuous(name=expression(beta),
                       breaks = 1:nrow(dat),labels = dat$beta)+
    scale_y_continuous(name='Power', #y轴坐标名称
                       breaks=seq(0,1,0.1),
                       limits=c(0,1))+
    theme(axis.title.x =element_text(size=14), 
          axis.title.y=element_text(size=14))+
    # scale_fill_manual(breaks = c("tpSSU", "tpSKAT", "SSU","SKAT"),
    #                   values=c("blue","red","green","violet"))+
    theme(legend.position="bottom")+#图例位置
    theme(legend.text = element_text(size = 10))+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))
  return(a)
}



