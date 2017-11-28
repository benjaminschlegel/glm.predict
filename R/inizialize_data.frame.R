initialize_data.frame = function(data, position = position){
  cnames = colnames(data)
  log.pos = grep("log\\(.+\\)",cnames)
  cnames = gsub("log\\(","",cnames)
  cnames = gsub("I\\(","",cnames)
  cnames = gsub("\\)","",cnames)
  cnames = gsub("\\^[0-9]+","",cnames)
  cnames = unique(cnames)
  if(!is.null(position)){
    dc_var = cnames[position]
    cnames[position] = paste0(dc_var,"_val1")
    cnames = append(cnames, paste0(dc_var,"_val2"), position)
    names = c("mean_val1","lower_val1","upper_val1","mean_val2","lower_val2","upper_val2","dc","dc_lower","dc_upper",cnames)
    log.pos = ifelse(log.pos > position, log.pos + 1, log.pos)
    if(position %in% log.pos){
      log.pos = append(log.pos, position + 1, which(log.pos == position))
    } 
    log.pos = log.pos + 9
  }else{
    names = c("mean","lower","upper",cnames)
    log.pos = log.pos + 3
  }
  result = data.frame(t(rep(NA,length(names))))
  colnames(result) = names
  return(list(result=result,log.pos=log.pos))
}
