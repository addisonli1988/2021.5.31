options(repos<- c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/") )
options("BioC_mirror"<- "https://mirrors.ustc.edu.cn/bioc/")

install.packages(c("devtools","curl")) ##Installs devtools and the MCPcounter dependancy 'curl'
library(devtools)
install_github("ebecht/MCPcounter",ref="master", subdir="Source")

library(MCPcounter)

estimate <- MCPcounter.estimate(rnaExpr, featuresType= "ENSEMBL_ID")
write.csv(estimate,"C:\\Users\\admin\\Desktop\\estimate.csv")
                                