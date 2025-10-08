#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")


#引用包
library(TwoSampleMR)

inputFile="GCST90277388.tsv.gz"      #输入文件
outFile="GCST90277388.exposure_data.csv"     #输出的结果文件
setwd("C:\\Users\\lexb\\Desktop\\lipidomeMedMR\\03.exposure")     #设置工作目录

#读取输入文件
rt=data.table::fread(inputFile, header=T, sep="\t", check.names=F)
rt$p_value=10^(-rt$neg_log_10_p_value)

#根据pvalue<1e-05对数据进行过滤(关联性分析)
data<-subset(rt, p_value<1e-05)
write.csv(file="exposure.pvalue.csv", data, row.names=F)

#读取输入文件, 并对输入文件进行格式转换
exposure_dat<-read_exposure_data(filename="exposure.pvalue.csv",
                       sep = ",",
                       snp_col = "rsid",
                       beta_col = "beta",
                       se_col = "standard_error",
                       effect_allele_col = "effect_allele",
                       other_allele_col = "other_allele",
                       eaf_col = "effect_allele_frequency",
                       pval_col = "p_value",
                       chr_col="chromosome",
                       pos_col = "base_pair_location",
                       clump = F)

#去除连锁不平衡的SNP
file.remove("exposure.pvalue.csv")
exposure_dat_clumped <- clump_data(exposure_dat, clump_kb=10000, clump_r2=0.001)
write.csv(exposure_dat_clumped, file=outFile, row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


