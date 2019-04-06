#R语言基础知识的一些检验，最好是对照几本R基础语法书籍来理解。

#全部答案及视频在：<https://github.com/jmzeng1314/R_bilibili>
  
#首先做完了周末班工作， 包括软件安装以及R包安装：<http://www.bio-info-trainee.com/3727.html>
  
#1.打开 `Rstudio` 告诉我它的工作目录。
getwd()

#2.新建6个向量，基于不同的`原子类型`。（重点是字符串，数值，逻辑值）
b<-"What's wrong"
aa<-1:10
c<-TRUE


#3.告诉我在你打开的rstudio里面 `getwd()` 代码运行后返回的是什么？

#当前工作目录

#4.新建一些数据结构，比如矩阵，数组，数据框，列表等(重点是数据框，矩阵）
newmatrix<-matrix(1:6, nrow = 2, ncol = 3)
newarray <- array(1:3, c(2,4))
pa1<-c(1,2,3,4)
ga1<-c(5,4,44,12)
pa2<-c(1,2,3,4)
ga2<-c(5,4,44,12)
pa3<-c(1,2,3,4)
ga3<-c(5,4,44,12)
newdata.frame <- data.frame(pa1,ga1,pa2,ga2,pa3,ga3)
newlist <- list(x = cars[,1], y = cars[,2])
newfactor <- factor(letters,ordered=TRUE)

                           
#5.在你新建的数据框进行切片操作，比如首先取第1，3行， 然后取第4，6列
newdata.frame[1:3,4:6]
                           
#6.使用data函数来加载R内置数据集 `rivers` 描述它。并且可以查看更多的R语言内置的数据集：<https://mp.weixin.qq.com/s/dZPbCXccTzuj0KkOL7R31g>
data("rivers")

                           
#7.下载 <https://www.ncbi.nlm.nih.gov/sra?term=SRP133642> 里面的 `RunInfo Table` 文件读入到R里面，了解这个数据框，多少列，每一列都是什么属性的元素。（参考B站生信小技巧获取runinfo table） 这是一个单细胞转录组项目的数据，共768个细胞，如果你找不到`RunInfo Table` 文件，可以[点击下载](http://www.bio-info-trainee.com/tmp/5years/SraRunTable.txt)，然后读入你的R里面也可以。

SraRunTable <- read.table("http://www.bio-info-trainee.com/tmp/5years/SraRunTable.txt",fill=TRUE,header = T,sep = "\t")
dim(SraRunTable)
class(colnames(SraRunTable))

                           
#8.下载 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111229> 里面的`样本信息sample.csv`读入到R里面，了解这个数据框，多少列，每一列都是什么属性的元素。（参考 <https://mp.weixin.qq.com/s/fbHMNXOdwiQX5BAlci8brA> 获取样本信息sample.csv）如果你实在是找不到样本信息文件sample.csv，也可以[点击下载](http://www.bio-info-trainee.com/tmp/5years/sample.csv)。
sample <- read.csv("sample.csv")

#9.把前面两个步骤的两个表（`RunInfo Table` 文件，样本信息sample.csv）关联起来，使用merge函数。
                           

d=merge(SraRunTable,sample,by.x = 'Sample_Name',by.y = 'Accession')# 在merge这里堵了很长时间...最终是数据源的问题 # by.x = 'Sample_Name',by.y = 'Accession'这里的Sample_Name与Accession必须代表同样的数据！
                           ```
### 基于下午的统计可视化
                           
#1.对前面读取的 `RunInfo Table` 文件在R里面探索其MBases列，包括 箱线图(boxplot)和五分位数(fivenum)，还有频数图(hist)，以及密度图(density) 。
                           

boxplot(MBases ~ MBytes, data = SraRunTable, col = "lightgray")
                           fivenum(SraRunTable$MBases)
                           plot(density(SraRunTable$MBases))
                           hist(SraRunTable$MBases)
                           plot(density(SraRunTable$MBases))
                           # ?fivenum
                           
#2.把前面读取的`样本信息`表格的样本名字`根据下划线分割`看第3列元素的统计情况。第三列代表该样本所在的plate
                           

e <- d[,c("MBases","Title")]
plate=unlist(lapply(as.character(e[,2]),function(x){# x=e[,2]  
x
strsplit(x,'_')[[1]][3]
}))
table(plate)
boxplot(e[,1]~plate)

                           
# 3.根据plate把关联到的 `RunInfo Table` 信息的MBases列分组检验是否有统计学显著的差异。
t.test(e[,1]~plate)

                           
#4.分组绘制箱线图(boxplot)，频数图(hist)，以及密度图(density) 。使用ggplot2把上面的图进行重新绘制。
                           
e$plate=plate
library(ggplot2)
colnames(e)
ggplot(e,aes(x=plate,y=MBases))+geom_boxplot()

                           
# 5.使用ggpubr把上面的图进行重新绘制。
                           

library(ggpubr)
p <- ggboxplot(e, x = "plate", y = "MBases",
               color = "plate", palette = "jco",
               add = "jitter")
               # Add p-value
p + stat_compare_means(method = 't.test')

                           
#6.随机取384个MBases信息，跟前面的两个plate的信息组合成新的数据框，第一列是分组，第二列是MBases，总共是384*3行数据。
                           

indexes <- sample (nrow (e), 384)
data <- e [indexes, ]
 data2 <- data[,c(3,1,2)]
dim(data2)
                           
#参考来源：生信技能树