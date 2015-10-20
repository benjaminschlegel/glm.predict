discrete.changes <-
function(model, values, sim.count=1000, conf.int=0.95){
  if(!is.character(values)){
    stop("values must be given as character!")
  }
  
  # check if any interaction
  formula = formula(model)
  dvs = unlist(strsplit(as.character(formula)[3]," + ", fixed=T))
  ia = c()
  for(dv in dvs){
    if(grepl(":",dv)){
      ia = c(ia,2)
      stop("Formula contains :, works only with * for interaction")
    }else if(grepl("\\*",dv)){
      ia = c(ia,1)
    }else{
      ia = c(ia, 0)
    }
  }
  
  formula = gsub(" \\* ", " + ", as.character(formula)[3])
  formula = gsub(" : ", " + ", formula)
  formula = gsub("\\*", " + ", formula)
  formula = gsub(":", " + ", formula)
  
  # remove dublicates
  temp.formula = unlist(strsplit(formula," + ", fixed=T))
  temp.formula = unique(temp.formula)
  
  var.pos = 1
  current.values = NA
  values.vector = unlist(strsplit(values,";"))
  k = 1
  global.counter = 1
  temp.results = list()
  factor.collector = data.frame(pos=NA,length=NA)
  factor.collector.counter = 1
  for(value in values.vector){
    #typen herausfinden
    multi = FALSE
    stopit = FALSE
    if(grepl("^mode$",value,ignore.case = TRUE)){
      temp = unlist(strsplit(as.character(formula)," \\+ "))[k]
      if(!is.null(model$data)){
        data = model$data
        variable = data[,grep(temp,colnames(data),value=T)[1]]
        temp2 = unique(variable)
        temp3 = temp2[which.max(table(variable))]
        if(is.numeric(temp3)){
          current.values = temp3
          stopit = TRUE
          global.counter = global.counter + 1
        }else{
          n = length(levels(variable))
          i = which(temp3==temp2)
          value = paste0("F",n,"(",i,")")
        }
        
      }
      
    }
    if(stopit){}
    else if(grepl("^F[0-9]+\\([0-9]+\\)$",value,ignore.case = TRUE)){
      temp = as.numeric(unlist(strsplit(value,"[F\\(\\)]")))
      if(temp[2]>2){
        multi = TRUE
      }
      n = temp[2]
      factor.collector[factor.collector.counter,]$pos = global.counter
      factor.collector[factor.collector.counter,]$length = n-1
      global.counter = global.counter + n-1
      factor.collector.counter = factor.collector.counter + 1
      for(i in 2:n){
        if(i==2){
          if(i == temp[3]){
            current.values = 1
          }else{
            current.values = 0
          }
        }else{
          if(i == temp[3]){
            current.values = c(current.values,1)
          }else{
            current.values = c(current.values,0)
          }
        }
      }
      current.values = t(t(current.values))
    }else if(grepl("^F[0-9]+$",value,ignore.case = TRUE)){
      temp = as.numeric(unlist(strsplit(value,"F")))[2]
      if(temp>2){
        multi = TRUE
      }
      n = temp[1]
      factor.collector[factor.collector.counter,]$pos = global.counter
      factor.collector[factor.collector.counter,]$length = n-1
      global.counter = global.counter + n-1
      factor.collector.counter = factor.collector.counter + 1
      
      high = temp-1
      current.values = matrix(nrow=high,ncol=2)
      for(j in 1:high){
        current.values[j,1] = 0
        current.values[j,2] = 1
      }
    }else if(grepl("^mean$",value,ignore.case = TRUE)){
      temp = unlist(strsplit(as.character(formula)," \\+ "))[k]
      if(!is.null(model$data)){
       data = model$data
       variable = data[,grep(temp,colnames(data),value=T)[1]]
       if(!is.numeric(variable)){
         stop("Cannot callculted the mean of a non numeric variable")
       }
       current.values = mean(variable, na.rm=T)
       global.counter = global.counter + 1
      }
    }else if(grepl("^median$",value,ignore.case = TRUE)){
      temp = unlist(strsplit(as.character(formula)," \\+ "))[k]
      if(!is.null(model$data)){
        data = model$data
        variable = data[,grep(temp,colnames(data),value=T)[1]]
        if(!is.numeric(variable)){
          stop("Cannot callculted the median of a non numeric variable")
        }
        current.values = median(variable, na.rm=T)
        global.counter = global.counter + 1
      }
    }else if(grepl("^[0-9\\.]+-[0-9\\.]+,[0-9[:punct:]]+$",value)){
      temp = as.numeric(unlist(strsplit(value,"[-,]")))
      current.values = seq(from=temp[1],to=temp[2],by=temp[3])
      global.counter = global.counter + 1
    }else if(grepl("^[0-9]+-[0-9]+$",value)){
      temp = unlist(strsplit(value,"-"))
      current.values = temp[1]:temp[2]
      global.counter = global.counter + 1
    }else if(grepl("^[0-9]+(,[0-9]+)+$",value)){
      current.values = as.numeric(unlist(strsplit(value,",")))
      global.counter = global.counter + 1
    }else if(grepl("^[0-9]+$",value)){
      current.values = as.numeric(value)
      global.counter = global.counter + 1
    }else{
      print(value)
      stop("values has invalid syntax!")
    }
    if(!multi){
      temp.results = c(temp.results, list(current.values))
    }else{
      c.values = list()
      for(i in 1:dim(current.values)[1]){
        c.values = c(c.values, list(current.values[i,]))
      }
      temp.results = c(temp.results, c.values)
      
    }
    #print(current.values)
    k = k + 1
  }
  
  temp = unlist(temp.results[var.pos])
  if(is.null(dim(temp))){
    n = length(temp)
  }else{
    n = dim(temp)[2]
  }
  
  product = 1
  products = c(1)
  i = 1
  in.factor = FALSE
  length = 0
  for(e in temp.results){
    if(in.factor && i > current.end){
      in.factor = FALSE
    }
    if(i != var.pos && !in.factor){
      if(is.null(dim(e))){
        n.temp = length(e)
      }else{
        n.temp = dim(e)[2]
      }
      products = c(products,product * n.temp)
      product = product * n.temp 
    }else if(length(which(factor.collector$pos==i))>0){
      row = which(factor.collector$pos==i)
      l = factor.collector[row,]$length
      length = l
      n.temp = 0
      for(k in 1:l){
        n.temp = n.temp + k
      }
      products = c(products,product * n.temp)
      product = product * n.temp
      in.factor = TRUE
      current.end = row + l - 1
    }
    i = i + 1
  }
  
  
  n.variables = length(temp.results)
  values1 = c()

  result = data.frame(mean1=NA,mean2=NA,lower1=NA,upper1=NA,lower2=NA,upper2=NA,mean.diff=NA,lower.diff=NA,upper.diff=NA)
  names = colnames(result)
  names.temp = unlist(strsplit(as.character(formula)," \\+ "))
  var.names = names.temp
  pos.found = FALSE
  for(i in 1:length(names.temp)){
    if(i == var.pos){
      old.name = names.temp[i]
      names.temp[i] = paste0(old.name,"1")
      if(!is.na(names.temp[i+1])){
        temp.save = names.temp[i+1]
      }
      names.temp[i+1] = paste0(old.name,"2")
      pos.found=TRUE
    }
    else if(pos.found){
      temp.old = temp.save
      temp.save = names.temp[i+1]
      names.temp[i+1] = temp.old
    }
  }
  result = cbind(result,data.frame(t(rep(NA,length(names.temp)))))
  names = c(names,names.temp)
  colnames(result) = names
  var.count = length(names.temp)
  
  factor.position = 0
  c.temp = 0
  c1 = 0
  c2 = 1
  c3 = 1
  count.checker = 0
  e.new = 0
  start_c3 = 1
  
  #print(products)
  
  for(i in 2:n){
    for(counter in 1:product){
      count.checker = count.checker + 1
      temp1 = c(1)
      temp2 = c(1)
      #print(temp.results)
      in.factor = FALSE
      colcounter = 10
      varcounter = 0
      did_varcount = FALSE
      t1_counter = 0
      t2_counter = 0
      do_colcounter = 0
      for(j in 1:n.variables){
        if(j != var.pos && !in.factor){
          if(j > var.pos){
            n.temp = length(unlist(temp.results[j]))
            prod_start = 1
            if(length(which(factor.collector$pos==j))>0){
              prod_start = j-1
            }
#             if(length(which(factor.collector$pos<j))>0){
#               rows = which(factor.collector$pos<j)
#               amount = sum(factor.collector[rows,]$length) - length(rows)
#               amount = 1
#             }else{
#               amount = 1
#             }
            e = (ceiling((counter/(products[j-1]/products[prod_start])))-1) %% n.temp + 1
          }else{
            n.temp = length(unlist(temp.results[j]))
            e = (ceiling((counter/products[j]))-1) %% n.temp + 1
          }
          t1 = unlist(temp.results[j])[e]
          t2 = unlist(temp.results[j])[e]
          if(length(which(factor.collector$pos<=j && factor.collector$pos+factor.collector$length>j))>0){
            row = which(factor.collector$pos<=j && factor.collector$pos+factor.collector$length>j)
            from = factor.collector[row,]$pos
            l = factor.collector[row,]$length
            to = from + l
            if(j == to-1){
              varcounter = varcounter + 1
              levels = levels(model$data[,var.names[varcounter]])
              if(sum(temp1[(from+1):(to-1)])+t1==0){
                result[counter+(i-2)*product,colcounter] = levels[1]
              }else if(t1 == 1){
                result[counter+(i-2)*product,colcounter] = levels[length(levels)]
              }else{
                for(p in 1:(length(levels)-2)){
                if(temp1[(from+p)]==1){
                    result[counter+(i-2)*product,colcounter] = levels[p+1]
                  }
                }
              }
              colcounter = colcounter + 1
            }
          }else{
            varcounter = varcounter + 1
            result[counter+(i-2)*product,colcounter] = t1
            colcounter = colcounter + 1
          }
          
        }else{
          if(length(which(factor.collector$pos<=j && factor.collector$pos+factor.collector$length>j))>0){
            row = which(factor.collector$pos<=j && factor.collector$pos+factor.collector$length>j)
            from = factor.collector[row,]$pos
            l = factor.collector[row,]$length
            to = from + l
            if(factor.position == 0){
              factor.position = j
            }
            
            if(count.checker>products[factor.position+1]){
              c.temp = 0
              c1 = 0
              c2 = 1
              c3 = 1
              start_c3 = 1
              e.new = 0
              count.checker = count.checker - products[factor.position+1]
            }
            
            e.old = e.new
            e.new = (ceiling((counter/products[factor.position]))-1) + 1
            
            c = floor(sqrt(products[factor.position+1] / products[factor.position] * 2))
            if(c.temp == 0){
              c.temp = c
            }
            
            if(c1<c){
              c1 = c1 + 1
            }else{
              c1 = 1
            }
            
            if(e.old != 0 && e.old != e.new){
              if(c2<c.temp){
                c2 = c2 + 1
              }else{
                c2 = 1
                c.temp = c.temp - 1
              }
              if(c3==c){
                start_c3 = start_c3 + 1
                c3 = start_c3
              }else{
                c3 = c3 + 1
              }
            }
            
            if(c-c.temp==c1){
              t1 = 1
              if(!did_varcount){
                varcounter = varcounter + 1
                did_varcount = TRUE
              }
              levels = levels(model$data[,var.names[varcounter]])
              result[counter+(i-2)*product,colcounter] = levels[c1+1]
              do_colcounter = do_colcounter + 1
              
            }else{
              t1 = 0
              t1_counter = t1_counter + 1
            }
            
            if(t1_counter == c){
              if(!did_varcount){
                varcounter = varcounter + 1
                did_varcount = TRUE
              }
              levels = levels(model$data[,var.names[varcounter]])
              result[counter+(i-2)*product,colcounter] = levels[1]
              do_colcounter = do_colcounter + 1
            }
            
            if(c3==c1){
              t2 = 1
              if(!did_varcount){
                varcounter = varcounter + 1
                did_varcount = TRUE
              }
              levels = levels(model$data[,var.names[varcounter]])
              result[counter+(i-2)*product,colcounter+1] = levels[c1+1]
              do_colcounter = do_colcounter + 1
            }else{
              t2 = 0
              t2_counter = t2_counter + 1
            }
            
            if(t2_counter == c){
              if(!did_varcount){
                varcounter = varcounter + 1
                did_varcount = TRUE
              }
              levels = levels(model$data[,var.names[varcounter]])
              result[counter+(i-2)*product,colcounter+1] = levels[1]
              do_colcounter = do_colcounter + 1
            }
            
            if(j < to-1){
              in.factor = TRUE
            }else{
              in.factor = FALSE
            }
            
            if(do_colcounter == 2){
              colcounter = colcounter + 2
              do_colcounter = 3
            }
          }else{
            t1 = unlist(temp.results[j])[i-1]
            t2 = unlist(temp.results[j])[i]
            result[counter+(i-2)*product,colcounter] = t1
            result[counter+(i-2)*product,colcounter+1] = t2
            colcounter = colcounter + 2
            varcounter = varcounter + 1
          }
          
        }
          temp1 = c(temp1,t1)
          temp2 = c(temp2,t2)
      }
      
      
      # test if non possible result
      c = 1
      ok = TRUE
      for(j in 1:var.count){
        if(length(which(factor.collector$pos==c))>0){
          row = which(factor.collector$pos==c)
          l = factor.collector[row,]$length
          from = c+1
          to = c+l
          if(sum(temp1[from:to])>1 | sum(temp2[from:to])>1){
            ok = FALSE
          }
          c = c + l
        }else{
          c = c + 1
        }
      }
      
      # interaction
      ia.counter = 2
      kick.out = c()
      for(ia.value in ia){
        if(ia.value==0){ #: +
          ia.counter = ia.counter + 1
        }
        if(ia.value==1 || ia.value==2){ # *
          if(ia.counter != 2){
            if(length(which(factor.collector$pos==ia.counter-1))>0){
              row = which(factor.collector$pos==ia.counter-1)
              l1 = factor.collector[row,]$length
            }else{
              l1 = 1
            }
          }else{
            l1 = 1
          }
          if(length(which(factor.collector$pos==ia.counter))>0){
            row = which(factor.collector$pos==ia.counter)
            l2 = factor.collector[row,]$length
          }else{
            l2 = 1
          }
          for(l1.i in 1:l1){
            for(l2.i in 1:l2){
              temp1 = c(temp1,temp1[ia.counter-1+l1.i]*temp1[ia.counter-1+l1+l2.i])
              temp2 = c(temp2,temp2[ia.counter-1+l1.i]*temp2[ia.counter-1+l1+l2.i])
            }
          }  
          if(ia.value==2){ # :
            kick.out = c(kick.out,(ia.counter):(ia.counter-1+l1+l2))
          }
        }
      }
      if(!is.null(kick.out)){
        temp1 = temp1[-(kick.out)]
      }
      
      if(ok){
         # print(temp1)
        #  print(temp2)
         # print("-------------------")
        subresult = discrete.change(model,temp1,temp2,sim.count,conf.int)
        result[counter+(i-2)*product,]$mean1 = subresult[1,1]
        result[counter+(i-2)*product,]$mean2 = subresult[2,1]
        result[counter+(i-2)*product,]$lower1 = subresult[1,2]
        result[counter+(i-2)*product,]$upper1 = subresult[1,3]
        result[counter+(i-2)*product,]$lower2 = subresult[2,2]
        result[counter+(i-2)*product,]$upper2 = subresult[2,3]
        result[counter+(i-2)*product,]$mean.diff = subresult[3,1]
        result[counter+(i-2)*product,]$lower.diff = subresult[3,2]
        result[counter+(i-2)*product,]$upper.diff = subresult[3,3]

      }
    }
  }
  #print(factor.collector)
  result = na.omit(result)
  rownames(result) = 1:length(result$mean1)
  return(result)
}
