pvalFilter=0.01      #pvalue过滤条件
mrFile="table.MRresult.csv"           #孟德尔随机化分析的结果文件
pleFile="table.pleiotropy.csv"        #多效性的结果文件
urlFile="lipidomeExposureID.txt"      #脂质体暴露数据ID的文件
setwd("C:\\Users\\lexb\\Desktop\\lipidomeMedMR\\07.lipidIVWfilter")     #设置工作目录

#读取孟德尔随机化的结果文件
rt=read.csv(mrFile, header=T, sep=",", check.names=F)

#根据pvalue对孟德尔随机化分析的结果进行过滤
ivwRT=rt[rt$method=="Inverse variance weighted",]
ivwRT=ivwRT[(ivwRT$pval<pvalFilter),]

#提取五种方法OR值方向一致的脂质体
ivw=data.frame()
for(lipidome in unique(ivwRT$exposure)){
	expoData=rt[rt$exposure==lipidome,]
	if(sum(expoData$or>1)==nrow(expoData) | sum(expoData$or<1)==nrow(expoData)){
		ivw=rbind(ivw, ivwRT[ivwRT$exposure==lipidome,])
	}
}

#读取多效性的结果文件
pleRT=read.csv(pleFile, header=T, sep=",", check.names=F)
#剔除pvalue小于0.05的脂质体
pleRT=pleRT[pleRT$pval>0.05,]
metLists=as.vector(pleRT$exposure)
outTab=ivw[ivw$exposure %in% metLists,]
row.names(outTab)=outTab[,1]

#在结果中添加脂质体数据的下载链接
urlRT=read.table(urlFile, header=T, sep="\t", check.names=F, row.names=1,  comment.char="", quote="")
sameID=intersect(row.names(outTab), row.names(urlRT))
outTab=cbind(outTab[sameID,,drop=F], urlRT[sameID,"URL",drop=F])
write.csv(outTab, file="lipidome-disease.IVWfilter.csv", row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


