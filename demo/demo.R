library(Butte)

load("exampleData.rda")

#an example with SCNA configuration 6:2 (total copy:minor copy)
test = Butte(x=exampleData_CN62$x, m=exampleData_CN62$m, history=cnmutHistory(6,2), nt=6, nb=2, qmethod="fullMLE",
      type = "butte", bootstrapCI="bootstrap", purity=exampleData_CN62$purity, B=100)
test$pi
test$piCI

#an example with SCNA configuration 4:1
test2 = Butte(x=exampleData_CN41$x, m=exampleData_CN41$m, history=cnmutHistory(4,1), nt=4, nb=1, qmethod="fullMLE",
      type = "identifiable", bootstrapCI="parametric", purity=exampleData_CN41$purity)
test2$pi
test2$piCI
