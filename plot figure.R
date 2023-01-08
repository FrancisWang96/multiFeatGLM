library(Metrics)

##########################置信区间作图######################################
xulie<-c(1:46)
#pre<-read.csv('/home/base/mrmr400_com_smile_xgb_lo_pre.csv',header=T)
range<-data.frame(sample=xulie,lab=pre$V4,igmu=pre$pre_mu,igup=pre$V3,igdown=pre$V2)
pic<-ggplot(data=range,aes(x=sample,y=lab))
pic<-pic+geom_ribbon(aes(ymin=igdown,ymax=igup,x=sample,fill="95% of ga"),alpha=0.4)
pic<-pic+geom_line(aes(x=xulie,y=igmu,colour="Prediction of ga"),size=1)
pic<-pic+geom_line(aes(x=xulie,y=igup),size=0.5,alpha=0.7)
pic<-pic+geom_line(aes(x=xulie,y=igdown),size=0.5,alpha=0.7)
pic<-pic+geom_hline(aes(yintercept=0),colour="#990000",linetype="dashed")
pic<-pic+scale_color_manual(values=c("blue4","deeppink4","lightskyblue","lightcoral"))
pic<-pic+geom_line(aes(x=xulie,y=lab,colour="True"),size=1)
pic<-pic+theme(panel.background = element_rect(fill="transparent",color="gray"))
pic<-pic+labs(title="prediction")
pic
#ggsave(pic,filename="mrmr400_com_smile_xgb_lo_pre.png",width=12,height=9)
############################################################################
#pre<-as.matrix(pre)
#pre[1,"V4"]

cp1<-cpi(pre1,220)
cp2<-cpi(pre2,46)
cp3<-cpi(pre3,42)
cp4<-cpi(pre4,45)
cp5<-cpi(pre5,46)
cp<-cpi(pre,220)
###########################置信区间覆盖率###################################
cpi<-function(pre,n){
  num=0
  for (i in 1:n) num=num+(pre[i,"V4"]>=pre[i,"V3"] & pre[i,"V4"]<=pre[i,"V2"])
  return(num/n)
}
############################################################################

###########################R Squared作图####################################

library(ggplot2)
lm_eqn<-function(pre){
  m<-lm(V4~pre_mu,pre);
  eq<-substitute(italic(r)^2~"="~r2,
                 list(r2=format(summary(m)$r.squared,digits=3)))
  as.character(as.expression(eq));
}
p<-ggplot(pre,aes(x=V4,y=pre_mu))+geom_point()
p+stat_smooth(method = 'lm',formula=y~x+I(x),colour='blue')+geom_text(x=-1,y=2,label=lm_eqn(pre),parse=TRUE)+theme(axis.test=element_text(colour='black',size=12),axis.title = element_text(size=14))+labs(x="Observed",y="Predicted")+theme_bw()+theme(panel.grid = element_blank())
p<-ggplot(pre,aes(x=V4,y=pre_mu))+geom_point()
p<-p1+stat_smooth(method = 'lm',formula=y~x,colour='blue')+labs(title="all prediction")+geom_text(x=10,y=20,label=lm_eqn(pre),parse=TRUE)+theme(axis.test=element_text(colour='black',size=12),axis.title = element_text(size=14))+labs(x="Observed",y="Predicted")+theme_bw()+theme(panel.grid = element_blank())

library(Rmisc)
xulie1<-c(1:41)
xulie2<-c(1:46)
xulie3<-c(1:42)
xulie4<-c(1:45)
xulie5<-c(1:46)
pre1<-read.csv('/home/gulab-01/妗岄潰/xgb_com_ga_epi_pre.csv',header=T)
pre2<-read.csv('/home/gulab-01/妗岄潰/xgb_com_ga_cis_pre.csv',header=T)
pre3<-read.csv('/home/gulab-01/妗岄潰/xgb_com_ga_cyc_pre.csv',header=T)
pre4<-read.csv('/home/gulab-01/妗岄潰/xgb_com_ga_doc_pre.csv',header=T)
pre5<-read.csv('/home/gulab-01/妗岄潰/xgb_com_ga_pac_pre.csv',header=T)

pre2<-rbind(pre1,pre2,pre3,pre4,pre5)
#pre2<-log(pre2)
rmse<-rmse(pre$pre_mu,pre$V4)
pre_mu<-log(pre$pre_mu)
log_y<-log(pre$V4)
rmse<-rmse(pre_mu,log_y)

pre_mu2<-log(pre2$pre_mu)
log_y2<-log(pre2$V4)
rmse2<-rmse(pre_mu2,log_y2)

pre_mu3<-log(pre3$pre_mu)
log_y3<-log(pre3$V4)
rmse3<-rmse(pre_mu3,log_y3)

pre_mu4<-log(pre4$pre_mu)
log_y4<-log(pre4$V4)
rmse4<-rmse(pre_mu4,log_y4)

pre_mu5<-log(pre5$pre_mu)
log_y5<-log(pre5$V4)
rmse5<-rmse(pre_mu5,log_y5)

p1<-ggplot(pre1,aes(x=V4,y=pre_mu))+geom_point()
p1<-p1+stat_smooth(method = 'lm',formula=y~x,colour='blue')+labs(title="epi prediction")+geom_text(x=10,y=20,label=lm_eqn(pre1),parse=TRUE)+theme(axis.test=element_text(colour='black',size=12),axis.title = element_text(size=14))+labs(x="Observed",y="Predicted")+theme_bw()+theme(panel.grid = element_blank())
#ggsave(p1,filename="mrmr400_com_smile_xgb_lo_rsquared.png",width=12,height=9)

p2<-ggplot(pre2,aes(x=V4,y=pre_mu))+geom_point()
p2<-p2+stat_smooth(method = 'lm',formula=y~x,colour='blue')+labs(title="cis prediction")+geom_text(x=505,y=1555,label=lm_eqn(pre2),parse=TRUE)+theme(axis.test=element_text(colour='black',size=12),axis.title = element_text(size=14))+labs(x="Observed",y="Predicted")+theme_bw()+theme(panel.grid = element_blank())

p3<-ggplot(pre3,aes(x=V4,y=pre_mu))+geom_point()
p3<-p3+stat_smooth(method = 'lm',formula=y~x,colour='blue')+labs(title="cyc prediction")+geom_text(x=1000,y=1000,label=lm_eqn(pre3),parse=TRUE)+theme(axis.test=element_text(colour='black',size=12),axis.title = element_text(size=14))+labs(x="Observed",y="Predicted")+theme_bw()+theme(panel.grid = element_blank())

p4<-ggplot(pre4,aes(x=V4,y=pre_mu))+geom_point()
p4<-p4+stat_smooth(method = 'lm',formula=y~x,colour='blue')+labs(title="dox prediction")+geom_text(x=5,y=5,label=lm_eqn(pre4),parse=TRUE)+theme(axis.test=element_text(colour='black',size=12),axis.title = element_text(size=14))+labs(x="Observed",y="Predicted")+theme_bw()+theme(panel.grid = element_blank())

p5<-ggplot(pre5,aes(x=V4,y=pre_mu))+geom_point()
p5<-p5+stat_smooth(method = 'lm',formula=y~x,colour='blue')+labs(title="pac prediction")+geom_text(x=2,y=40,label=lm_eqn(pre5),parse=TRUE)+theme(axis.test=element_text(colour='black',size=12),axis.title = element_text(size=14))+labs(x="Observed",y="Predicted")+theme_bw()+theme(panel.grid = element_blank())
pic<-grid.arrange(p1,p2,p3,p4,p5,p,ncol=2)
ggsave(pic,filename="com_xgb_ga_rsquare.png",width=12,height=9)
############################################################################

##########################置信区间作图######################################

range1<-data.frame(sample=xulie1,lab=pre1$V4,igmu=pre1$pre_mu,igup=pre1$V3,igdown=pre1$V2)
range2<-data.frame(sample=xulie2,lab=pre2$V4,igmu=pre2$pre_mu,igup=pre2$V3,igdown=pre2$V2)
range3<-data.frame(sample=xulie3,lab=pre3$V4,igmu=pre3$pre_mu,igup=pre3$V3,igdown=pre3$V2)
range4<-data.frame(sample=xulie4,lab=pre4$V4,igmu=pre4$pre_mu,igup=pre4$V3,igdown=pre4$V2)
range5<-data.frame(sample=xulie5,lab=pre5$V4,igmu=pre5$pre_mu,igup=pre5$V3,igdown=pre5$V2)

pic1<-ggplot(data=range1,aes(x=sample,y=lab))
pic1<-pic1+geom_ribbon(aes(ymin=igdown,ymax=igup,x=sample,fill="95% of ga"),alpha=0.4)
pic1<-pic1+geom_line(aes(x=xulie1,y=igmu,colour="Prediction of ga"),size=1)
pic1<-pic1+geom_line(aes(x=xulie1,y=igup),size=0.5,alpha=0.7)
pic1<-pic1+geom_line(aes(x=xulie1,y=igdown),size=0.5,alpha=0.7)
pic1<-pic1+geom_hline(aes(yintercept=0),colour="#990000",linetype="dashed")
pic1<-pic1+scale_color_manual(values=c("blue4","deeppink4","lightskyblue","lightcoral"))
pic1<-pic1+geom_line(aes(x=xulie1,y=lab,colour="True"),size=1)
pic1<-pic1+theme(panel.background = element_rect(fill="transparent",color="gray"))
pic1<-pic1+labs(title="epi prediction")

pic2<-ggplot(data=range2,aes(x=sample,y=lab))
pic2<-pic2+geom_ribbon(aes(ymin=igdown,ymax=igup,x=sample,fill="95% of ga"),alpha=0.4)
pic2<-pic2+geom_line(aes(x=xulie2,y=igmu,colour="Prediction of ga"),size=1)
pic2<-pic2+geom_line(aes(x=xulie2,y=igup),size=0.5,alpha=0.7)
pic2<-pic2+geom_line(aes(x=xulie2,y=igdown),size=0.5,alpha=0.7)
pic2<-pic2+geom_hline(aes(yintercept=0),colour="#990000",linetype="dashed")
pic2<-pic2+scale_color_manual(values=c("blue4","deeppink4","lightskyblue","lightcoral"))
pic2<-pic2+geom_line(aes(x=xulie2,y=lab,colour="True"),size=1)
pic2<-pic2+theme(panel.background = element_rect(fill="transparent",color="gray"))
pic2<-pic2+labs(title="cis prediction")

pic3<-ggplot(data=range3,aes(x=sample,y=lab))
pic3<-pic3+geom_ribbon(aes(ymin=igdown,ymax=igup,x=sample,fill="95% of ga"),alpha=0.4)
pic3<-pic3+geom_line(aes(x=xulie3,y=igmu,colour="Prediction of ga"),size=1)
pic3<-pic3+geom_line(aes(x=xulie3,y=igup),size=0.5,alpha=0.7)
pic3<-pic3+geom_line(aes(x=xulie3,y=igdown),size=0.5,alpha=0.7)
pic3<-pic3+geom_hline(aes(yintercept=0),colour="#990000",linetype="dashed")
pic3<-pic3+scale_color_manual(values=c("blue4","deeppink4","lightskyblue","lightcoral"))
pic3<-pic3+geom_line(aes(x=xulie3,y=lab,colour="True"),size=1)
pic3<-pic3+theme(panel.background = element_rect(fill="transparent",color="gray"))
pic3<-pic3+labs(title="cyc prediction")

pic4<-ggplot(data=range4,aes(x=sample,y=lab))
pic4<-pic4+geom_ribbon(aes(ymin=igdown,ymax=igup,x=sample,fill="95% of ga"),alpha=0.4)
pic4<-pic4+geom_line(aes(x=xulie4,y=igmu,colour="Prediction of ga"),size=1)
pic4<-pic4+geom_line(aes(x=xulie4,y=igup),size=0.5,alpha=0.7)
pic4<-pic4+geom_line(aes(x=xulie4,y=igdown),size=0.5,alpha=0.7)
pic4<-pic4+geom_hline(aes(yintercept=0),colour="#990000",linetype="dashed")
pic4<-pic4+scale_color_manual(values=c("blue4","deeppink4","lightskyblue","lightcoral"))
pic4<-pic4+geom_line(aes(x=xulie4,y=lab,colour="True"),size=1)
pic4<-pic4+theme(panel.background = element_rect(fill="transparent",color="gray"))
pic4<-pic4+labs(title="doc prediction")

pic5<-ggplot(data=range5,aes(x=sample,y=lab))
pic5<-pic5+geom_ribbon(aes(ymin=igdown,ymax=igup,x=sample,fill="95% of ga"),alpha=0.4)
pic5<-pic5+geom_line(aes(x=xulie5,y=igmu,colour="Prediction of ga"),size=1)
pic5<-pic5+geom_line(aes(x=xulie5,y=igup),size=0.5,alpha=0.7)
pic5<-pic5+geom_line(aes(x=xulie5,y=igdown),size=0.5,alpha=0.7)
pic5<-pic5+geom_hline(aes(yintercept=0),colour="#990000",linetype="dashed")
pic5<-pic5+scale_color_manual(values=c("blue4","deeppink4","lightskyblue","lightcoral"))
pic5<-pic5+geom_line(aes(x=xulie5,y=lab,colour="True"),size=1)
pic5<-pic5+theme(panel.background = element_rect(fill="transparent",color="gray"))
pic5<-pic5+labs(title="pac prediction")

library(gridExtra)
pic<-grid.arrange(pic1,pic2,pic3,pic4,pic5,pic,ncol=2)
ggsave(pic,filename="mrmr_com_xgb_ga_pre.png",width=12,height=9)
###############################################################################################