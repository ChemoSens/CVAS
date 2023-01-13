library(CVAS)
library()
data(cheeses)
cheesesR1=cheeses[cheeses[,"Replicate"]=="R1",]
tableProduct=aggregate(cheeses[,-c(1:4)],by=list(cheeses[,"ProductCode"]),FUN=mean)
attributes=colnames(cheeses)[-c(1:4)]

rescva_ow=CVA(cheeses,representation="twoMaps",option="ow",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
rescva_tw=CVA(cheeses,representation="twoMaps",option="tw",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
rescva_tws=CVA(cheeses,representation="twoMaps",option="tws",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
rescva_mam=CVA(cheeses,representation="twoMaps",option="mam",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
rescva_overall=CVA(cheeses,representation="twoMaps",option="overall",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
rescva_1w=CVA(cheeses,representation="twoMaps",option="1w",nbAxes=2,productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")

rescva_owb=CVA(cheeses,option="ow",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
rescva_twb=CVA(cheeses,option="tw",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
rescva_twsb=CVA(cheeses,option="tws",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
rescva_mamb=CVA(cheeses,option="mam",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
rescva_overallb=CVA(cheeses,option="overall",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
rescva_1wb=CVA(cheeses,option="1w",nbAxes=2,productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")

resManova=manova(as.matrix(cheeses[,-c(1:4)])~cheeses[,"ProductCode"])
summary(resManova,test="Wilks")
rescva_ow$Stats

resManova=manova(as.matrix(cheeses[,-c(1:4)])~cheeses[,"ProductCode"]+cheeses[,"SubjectCode"])
summary(resManova,test="Wilks")
rescva_tws$Stats

rescva_ow$B


# Comparison between cva 1w et lda
#======================================
reslda=lda(cheeses[,-c(1:4)],grouping=cheeses[,"ProductCode"])
test_that("Correlation of the weights of 1 between lda and cva 1w (axis1)",
          expect_true(abs(cor(reslda$scaling,rescva_1w$EigenVectors[,1:2])[1,1])==1)
          )
reslda=lda(cheeses[,-c(1:4)],grouping=cheeses[,"ProductCode"])
test_that("Correlation of the weights of 1 between lda and cva 1w (axis2)",
          expect_true(abs(cor(reslda$scaling,rescva_1w$EigenVectors[,1:2])[2,2])==1)
)

#/!\abs(cor(reslda$scaling,rescva_1wb$EigenVectors[,1:2]))[1,1]

# Comparison between cva tw and LDA on averaged means
#=====================================================
# cheesesAvg=aggregate(cheeses[,attributes],by=list(cheeses[,"ProductCode"],cheeses[,"SubjectCode"]),FUN="mean")
# colnames(cheesesAvg)[1:2]=c("Product","Subject")
# reslda=lda(cheesesAvg[,-c(1:2)],grouping=cheesesAvg[,"Product"])



# Checking biplot results (relation with tm)
 #==================================================
 test_that("Correlations between biplot results and two maps results (first axis)",
           expect_true(
             cor(rescva_tw$IndivCoord,rescva_tw$IndivCoord[,1:2])[1,1]==1
           ))
 test_that("Correlations between biplot results and two maps results (second axis)",
           expect_true(
             cor(rescva_tw$IndivCoord,rescva_twb$IndivCoord[,1:2])[2,2]==1
           ))
