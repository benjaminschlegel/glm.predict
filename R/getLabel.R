getLabel = function(data,i,pos){
  labels = levels(data[,i+1]) # +1 because of y
  return(labels[pos])
}