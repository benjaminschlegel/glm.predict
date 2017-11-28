getNames = function(names){
  new.names = c("mean","lower","upper")
  for(i in 1:length(names)){
    new.names = c(new.names,names[i])
  }
  result = data.frame(t(rep(NA,length(new.names))))
  colnames(result) = new.names
  return(result)
}