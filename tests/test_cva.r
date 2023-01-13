library(CVAS)
data(cheeses)
resCVA=CVA(cheeses,option="tw",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")
rescva1=CVA(cheeses,option="ow",productName="ProductCode",subjectName="SubjectCode",replicateName="Replicate",sessionName="Session")

#
#
#
# # Recovering same results as lda
# extendedDataRep1=extendedData[extendedData[,"rep"]==1,]
# rescva1=CVA(extendedDataRep1,option="OneWayANOVA")
# library(MASS)
# reslda=lda(extendedDataRep1[,-c(1:3)],grouping=extendedDataRep1[,"product"])
# test_that("Correlation of the weights of 1 between lda and cva (axis1)",
#           expect_true(
#             abs(cor(reslda$scaling,rescva1$EigenVectors[,1:2]))[1,1]==1
#           ))
#
# test_that("Correlation of the weights of 1 between lda and cva (axis2)",
#           expect_true(
#             abs(cor(reslda$scaling,rescva1$EigenVectors[,1:2]))[2,2]==1
#           ))
#
# predTest<-predict(reslda)
# predTest$x
# rescva1$IndivCoord
# rescva1_gg=turnToCVAgg(rescva1)
# test_that("Correlation 1 between individual scores of lda and cva (axis1)",
#           expect_true(
# abs(cor(predTest$x[,1],rescva1_gg$indSup[,1]))==1
# ))
#
# test_that("Correlation 1 between individual scores of lda and cva (axis 2)",
#           expect_true(
# abs(cor(predTest$x[,2],rescva1_gg$indSup[,2]))==1
#           ))
#
#
