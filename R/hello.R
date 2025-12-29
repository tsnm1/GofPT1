# install.packages("usethis", "devtools", "roxygen2")

library(usethis)
library(devtools)
library(roxygen2)
library(renv)

# devtools::has_devel()

# usethis::create_package(path = "E:/一些论文/假设检验/Knockoff方式/R_package/ZIPGSK")

# 初始化环境
renv::init()
renv::snapshot()

# 如果包为使用时必须的，则需要设置 type = "Imports"
# {renv} 在这里只是开发必备，而非使用所开发的包必备，因此选择 "Suggests" 即可
usethis::use_package(package = "renv", type = "Suggests")
# usethis::use_package(package = "renv", type = "Suggests")

# usethis::use_package(package = "knockoff", type = "Imports")
# usethis::use_package(package = "scDesign2", type = "Imports")
# usethis::use_package(package = "ZIPG", type = "Imports")

usethis::use_package(package = "MASS", type = "Imports")
usethis::use_package(package = "glmnet", type = "Imports")
usethis::use_package(package = "LassoSIR", type = "Imports")
usethis::use_package(package = "VariableScreening", type = "Imports")
usethis::use_package(package = "psych", type = "Imports")
usethis::use_package(package = "caret", type = "Imports")


# 引入必备包 {rmarkdown}
# install.packages("rmarkdown")
library(rmarkdown)

# 虽然也可以使用 usethis::use_readme_md()，看个人需求（不详细讲解区别）
usethis::use_readme_rmd()

# 写完后，需要生成文档才能够真正使用
devtools::document()

devtools::install()

# 安装 {styler}
# install.packages("styler")
library(styler)
# 对整个包进行代码美化
styler::style_pkg()


# data

data_crime <-  read.csv("D:/坚果云/一些论文/假设检验/Gof_test/数据集/communities+and+crime/communities.data", header = FALSE)
data_crime <- data_crime[, -c(1:5)]
data_crime[data_crime == "?"] <- NA
data_crime <- apply(data_crime, 2, as.numeric)
data_crime <- data_crime[, !is.na(apply(data_crime, 2, sum))]
# na_prop <- colMeans(is.na(data_crime))
# data_crime <- data_crime[, na_prop <= 0.5]
dim(data_crime)
sum(is.na(data_crime))
y <- data_crime[, 100]
x <- data_crime[, -100]
colnames(data_crime) <- c(paste("X",1:99,sep = ""),"Y")


data_AML <- read.csv("D:/坚果云/一些论文/假设检验/Gof_test/数据集/aml_ohsu_2022/data_AML.csv")
dim(data_AML)
y <- data_AML$ELN_binary
x <- data_AML[,1:(ncol(data_AML)-2)]
x <- x[,-which(is.na(apply(data_AML,2,sum)))]
# x <- apply(x,2,as.numeric)
sum(is.na(x))

data_full_filtered$`De Novo` <- ifelse(data_full_filtered$`De Novo` =="TRUE",1,0 )

usethis::use_data(data_crime,data_AML, overwrite = TRUE)

# getwd()
# load("E:/一些论文/假设检验/Knockoff方式/R_package/SKMSD/data/data_pg_copula.rda")


# 3.9 编写包的说明
usethis::use_mit_license()

# 3.10 用 {pkgdown} 制作包的说明书
# install.packages("pkgdown")
library(pkgdown)
# 初始化你的用户手册网站
usethis::use_pkgdown()

pkgdown::build_site()


# 3.11 用 Github Action 自动检查
usethis::use_git()
usethis::use_github()
usethis::use_tidy_github_actions()


# 安装包
# devtools::install_github("tsnm1/SKMSD")
# library(SKMSD)
#
# SKMSD
# data("data_pg_copula")
