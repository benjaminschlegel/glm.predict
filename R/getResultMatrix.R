getResultMatrix = function(part_y, levels, result_matrix, base.combinations){
  return(unlist(lapply(1:length(base.combinations[,1]), getMatrixPart, part_y = part_y, levels = levels, result_matrix = result_matrix)))
}

getMatrixPart = function(part_x, part_y, levels, result_matrix){
  from = levels * (part_y - 1) + 1
  to = from + levels - 1
  return(result_matrix[from:to, part_x])
}

