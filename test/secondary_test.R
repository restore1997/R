
  
### 作业 1
  
#请根据R包`org.Hs.eg.db`找到下面ensembl 基因ID 对应的基因名(symbol)
```
ENSG00000000003.13
ENSG00000000005.5
ENSG00000000419.11
ENSG00000000457.12
ENSG00000000460.15
ENSG00000000938.11
```
# 先将待解的ensembl基因ID生成txt，为e1.txt
rm(list = ls())
options(stringsAsFactors = F)
a=read.table('e1.txt')
head(a)
# 载入注释包（需提前安装，见下）
library(org.Hs.eg.db)
# 查看包内容
ls("package:org.Hs.eg.db")
g2s=toTable(org.Hs.egSYMBOL);head(g2s)
g2e=toTable(org.Hs.egENSEMBL);head(g2e)
head(g2e)
library(stringr)
# 用strsplit对a进行分割，取"."之前的内容，即ensembl_id
a$ensembl_id=unlist(lapply(a$V1,function(x){
  strsplit(as.character(x),'[.]')[[1]][1]
})
)
# 将三者关联
tmp1=merge(a,g2e,by='ensembl_id')
tmp2=merge(tmp1,g2s,by='gene_id')

### 作业 2
#根据R包`hgu133a.db`找到下面探针对应的基因名(symbol)
1053_at
117_at
121_at
1255_g_at
1316_at
1320_at
1405_i_at
1431_at
1438_at
1487_at
1494_f_at
1598_g_at
160020_at
1729_at
177_at

rm(list = ls())
options(stringsAsFactors = F)
# 将探针保存为txt，命名为e2，并读取
a=read.table('e2.txt')
# 修改列名
colnames(a)='probe_id'
library(hgu133a.db)
# 找到注释文件内探针与基因的对应信息
ids=toTable(hgu133aSYMBOL)
head(ids)
# 合并
tmp1=merge(ids,a,by='probe_id')
# tmp1和tmp2效果相同
tmp2=ids[match(a$probe_id,ids$probe_id),]
# 如果之前不修改a的列名，可用
tmp3=merge(ids,a,by.x='probe_id',by.y="V1")

### 作业 3
#找到R包`CLL`内置的数据集的表达矩阵里面的TP53基因的表达量，并且绘制在 `progres.-stable`分组的boxplot图。想想如何通过 `ggpubr` 进行美化。

rm(list = ls())
options(stringsAsFactors = F)
# 不提示加载信息
suppressPackageStartupMessages(library(CLL))
# 载入CLL信息
data(sCLLex)
sCLLex
# 获取表达矩阵及样本信息
exprSet=exprs(sCLLex)
pd=pData(sCLLex)
# 加载注释包
library(hgu95av2.db) 
ids=toTable(hgu95av2SYMBOL)
head(ids)
# 在ids直接搜索TP53，得三个探针；或者用代码获取（探针数量多时推荐）
aaa <- ids[ids$symbol%in%"TP53",][,1]
#画图
boxplot(exprSet['1939_at',] ~ pd$Disease)
boxplot(exprSet['1974_s_at',] ~ pd$Disease)
boxplot(exprSet['31618_at',] ~ pd$Disease)

### 作业6
#下载数据集GSE17215的表达矩阵并且提取下面的基因画热图

ACTR3B ANLN BAG1 BCL2 BIRC5 BLVRA CCNB1 CCNE1 CDC20 CDC6 CDCA1 CDH3 CENPF CEP55 CXXC5 EGFR ERBB2 ESR1 EXO1 FGFR4 FOXA1 FOXC1 GPR160 GRB7 KIF2C KNTC2 KRT14 KRT17 KRT5 MAPT MDM2 MELK MIA MKI67 MLPH MMP11 MYBL2 MYC NAT1 ORC6L PGR PHGDH PTTG1 RRM2 SFRP1 SLC39A6 TMEM45B TYMS UBE2C UBE2T

#提示：根据基因名拿到探针ID，缩小表达矩阵绘制热图，没有检查到的基因直接忽略即可。

rm(list = ls()) 
options(stringsAsFactors = F)
# 注意查看下载文件的大小，检查数据 
f='GSE17215_eSet.Rdata'
# 使用GEOquery包下载GEO数据（需提前下载安装此包）
library(GEOquery)
# 下载可能会出错，可以试试更改镜像、翻墙，重启，rm（list=ls（），据说早上的网好容易下载...看运气了，下2为更改镜像
options(repos<- c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/") )
options("BioC_mirror"<- "https://mirrors.ustc.edu.cn/bioc/")
if(!file.exists(f)){
  gset <- getGEO('GSE17215', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   ## 保存到本地
}
# 先保存后载入，养成好习惯
load('GSE17215_eSet.Rdata')
# 因为这个GEO数据集只有一个GPL平台，所以下载到的是一个含有一个元素的list
class(gset)
length(gset)
class(gset[[1]])
# not list易于操作
a=gset[[1]]
# 获取表达矩阵
dat=exprs(a)
dim(dat)
# 加载注释包
library(hgu133a.db)
# 看一下包里有什么
ls("package:hgu133a.db")
# 获取探针与基因的对应关系
ids=toTable(hgu133aSYMBOL)
head(ids)
# 筛选ids中包含的probe_id，将两个表行名设为同序
dat=dat[ids$probe_id,]
dat[1:4,1:4] 
# 取dat表达量的中位数，生成ids的median列
ids$median=apply(dat,1,median)
# 将ids按照symbol和median排序（降序），可以看一下前后对比关系进行理解
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
# 去除重复的symbol值，保留第一位，即表达量最大的一位（dim一下看看前后对比）
ids=ids[!duplicated(ids$symbol),]
# 同上
dat=dat[ids$probe_id,]
# 先同序，再直接改为symbol（dim一下瞧瞧两者行数）
rownames(dat)=ids$symbol
dat[1:4,1:4]  
dim(dat)
# 好，进入正题
ng='ACTR3B ANLN BAG1 BCL2 BIRC5 BLVRA CCNB1 CCNE1 CDC20 CDC6 CDCA1 CDH3 CENPF CEP55 CXXC5 EGFR ERBB2 ESR1 EXO1 FGFR4 FOXA1 FOXC1 GPR160 GRB7 KIF2C KNTC2 KRT14 KRT17 KRT5 MAPT MDM2 MELK MIA MKI67 MLPH MMP11 MYBL2 MYC NAT1 ORC6L PGR PHGDH PTTG1 RRM2 SFRP1 SLC39A6 TMEM45B TYMS UBE2C UBE2T'
# 把基因加上“”，使之与rownames(dat)格式一致
ng=strsplit(ng,' ')[[1]]
# 看一下，多少基因存在于dat中，可得41个TRUE&9个FALSE
table(ng %in%  rownames(dat))
# 清洗掉不存在的ng，注意这一步存在排序（连同下一步理解）
ng=ng[ng %in%  rownames(dat)]
# 再改行名
dat=dat[ng,]
# 画图
dat=log2(dat)
pheatmap::pheatmap(dat,scale = 'row')

### 作业7
#下载数据集GSE24673的表达矩阵计算样本的相关性并且绘制热图，需要标记上样本分组信息

rm(list = ls())  
options(stringsAsFactors = F)
GSE_name = 'GSE24673'
options( 'download.file.method.GEOquery' = 'libcurl' ) #windows系统
gset <- getGEO( GSE_name, getGPL = F )
# 保存并载入数据
save( gset, file = 'GSE24673_gset.Rdata' )
load(file ="GSE24673_gset.Rdata")
a <- gset[[1]]
dat=exprs(a)
dim(dat)
pd=pData(a)
# 我们通过观察样本信息（pd）来得知11个样本的性质（分组信息），建立小分组
group_list=c('rbc','rbc','rbc',
             'rbn','rbn','rbn',
             'rbc','rbc','rbc',
             'normal','normal')
dat[1:4,1:4]
# 相关性分析并作图
M=cor(dat)
pheatmap::pheatmap(M)
# 标记分组信息，相关性分析并作图
tmp=data.frame(g=group_list)
rownames(tmp)=colnames(M)
pheatmap::pheatmap(M,annotation_col = tmp)

### 作业8
#找到 GPL6244 platform of Affymetrix Human Gene 1.0 ST Array 对应的R的bioconductor注释包，并且安装它！

#去菜鸟团网站搜索 http://www.bio-info-trainee.com/1399.html
options()$repos
options()$BioC_mirror 
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
BiocManager::install("hugene10sttranscriptcluster.db",ask = F,update = F)
options()$repos
options()$BioC_mirror

### 作业9
#下载数据集GSE42872的表达矩阵，并且分别挑选出 所有样本的(平均表达量/sd/mad/)最大的探针，并且找到它们对应的基因。

# 经典步骤前面已经解释过很多次了
rm(list = ls())  
options(stringsAsFactors = F)
f='GSE42872_eSet.Rdata'
library(GEOquery)
if(!file.exists(f)){
  gset <- getGEO('GSE42872', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f) 
}
load('GSE42872_eSet.Rdata') 
class(gset)
length(gset)
class(gset[[1]])
a=gset[[1]]
dat=exprs(a)
dim(dat)
pd=pData(a)
# 先画个箱图看看
boxplot(dat)
# 分别挑选出 所有样本的(平均表达量/sd/mad/)最大的探针，分别为7978905，8133876和8133876
sort(apply(dat,1,mean),decreasing = T)[1]
sort(apply(dat,1,sd),decreasing = T)[1]
sort(apply(dat,1,mad),decreasing = T)[1]
# 加载bioconduct注释包（包名有没有很熟悉？）
library(hugene10sttranscriptcluster.db)
ls("package:hugene10sttranscriptcluster.db")
ggg <-toTable(hugene10sttranscriptclusterSYMBOL)
# 找一下；尴尬的是7978905找不到...
ggg2 <- ggg[ggg$probe_id%in%7978905,]
ggg2 <- ggg[ggg$probe_id%in%8133876,]

### 作业10
#下载数据集GSE42872的表达矩阵，并且根据分组使用limma做差异分析，得到差异结果矩阵

rm(list = ls()) 
options(stringsAsFactors = F)
f='GSE42872_eSet.Rdata'
library(GEOquery)
if(!file.exists(f)){
  gset <- getGEO('GSE42872', destdir=".",
                 AnnotGPL = F,     ## 注释文件
                 getGPL = F)       ## 平台文件
  save(gset,file=f)   
}
load('GSE42872_eSet.Rdata') 
class(gset)
length(gset)
class(gset[[1]])
a=gset[[1]]
dat=exprs(a)
dim(dat)
pd=pData(a)
boxplot(dat)
# 从样本信息提取出分组
group_list=unlist(lapply(pd$title,function(x){
  strsplit(x,' ')[[1]][4]
}))
exprSet=dat
exprSet[1:4,1:4]
# 用limma包做差异表达分析
#差异分析重点在于做好表达矩阵和分组信息，具体原理可以不用理解
suppressMessages(library(limma)) 
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix 
##step1
fit <- lmFit(exprSet,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)