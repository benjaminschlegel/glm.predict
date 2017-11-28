getLabel = function(data,i,pos){
  labels = levels(data[,i])
  return(labels[pos])
}