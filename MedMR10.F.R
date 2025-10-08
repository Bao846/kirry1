#install.packages("devtools")
#devtools::install_github("mrcieu/ieugwasr", force = TRUE)


#引用包
library(ieugwasr)

inputFile="met.exposure_data.csv"      #输入文件
setwd("C:\\Users\\lexb\\Desktop\\lipidomeMedMR\\10.metF")     #设置工作目录

#读取输入文件
dat=read.csv(inputFile, header=T, sep=",", check.names=F)

#计算F检验值
dat$R2<-(2*dat$beta.exposure*dat$beta.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)/(2*dat$beta.exposure*dat$beta.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)+2*dat$se.exposure*dat$se.exposure*dat$samplesize.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)))     #计算R2
dat$F<-dat$R2*(dat$samplesize.exposure-2)/(1-dat$R2)     #计算F检验值

#根据F值>10对数据进行过滤, 删除弱工具变量
outTab=dat[as.numeric(dat$F)>10,]
write.csv(outTab, file="met.exposure.F.csv", row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio

