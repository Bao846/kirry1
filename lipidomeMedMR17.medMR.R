setwd("C:\\Users\\lexb\\Desktop\\lipidomeMedMR\\17.medMR")      #设置工作目录
files=dir()                           #获取目录下所有文件
files=grep("csv$", files, value=T)    #提取csv结尾的文件
beta1File=grep("^beta1.", files, value=T)          #脂质体到代谢物MR分析结果文件
beta2File=grep("^beta2.", files, value=T)          #代谢物到疾病MR分析结果文件
betaAllFile=grep("^betaAll.", files, value=T)      #脂质体到疾病MR分析结果文件

#读取三个输入文件, 获取beta值
beta1RT=read.csv(beta1File, header=T, sep=",", check.names=F)
beta2RT=read.csv(beta2File, header=T, sep=",", check.names=F)
betaAllRT=read.csv(betaAllFile, header=T, sep=",", check.names=F)
b1=beta1RT[3,"b"]
se1=beta1RT[3,"se"]
b2=beta2RT[3,"b"]
se2=beta1RT[3,"se"]
bataAll=betaAllRT[3,"b"]

#计算中介效应(indirect effect)
beta12=b1 * b2
se12=sqrt((b1^2)*(se1^2)+(b2^2)*(se2^2))
ciLow=beta12-1.96*se12
ciHigh=beta12+1.96*se12
medEffect=paste0(signif(beta12, digits=3), "(", signif(ciLow, digits=3), ", ", signif(ciHigh, digits=3), ")")

#计算直接效应(direct effect)
betaDirect=bataAll-beta12
betaDirect

#中介效应所占的百分比
beta12Ratio=beta12/bataAll
ciLowRatio=ciLow/bataAll
ciHighRatio=ciHigh/bataAll
medRatioEffect=paste0(signif(beta12Ratio*100, digits=3), "%(", signif(ciLowRatio*100, digits=3), "%, ", signif(ciHighRatio*100, digits=3), "%)")
#medRatioEffect=paste0(signif(beta12Ratio*100, digits=3), "%")

#计算pvalue
zValue=beta12/se12
pvalue=2*pnorm(q=abs(zValue), lower.tail=FALSE)

#输出结果的表格
outTab=data.frame(unique(beta1RT[,"exposure"]),
				  unique(beta2RT[,"exposure"]),
				  unique(betaAllRT[,"outcome"]),
                  medEffect, medRatioEffect, pvalue)
colnames(outTab)=c("Lipidome","Metabolite","outcome","Mediated effect", "Mediated proportion", "pvalue")
write.csv(outTab, file="MediatedEffect.csv", row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


