
setwd("D:/Dropbox/VM 1")
setwd("D:/Dropbox/VM 2")
setwd("D:/UW Coursework/AWS Final Parameter Set 3600 results")

length(grep(list.files(),pattern="param",value=T))

merge_list<-paste0("parameter_set_",1,"_")
for (i in 2:64){merge_list<-c(merge_list,paste0("parameter_set_",i,"_"))}


for (i in merge_list){
  write.csv(Reduce(rbind,lapply(grep(list.files(),pattern=i,value=T),read.csv,header=F)),file=paste0(i,"final.csv"),row.names = F)
}